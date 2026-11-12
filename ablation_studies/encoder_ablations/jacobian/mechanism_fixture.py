"""Tiny CPU-only native-network mechanism checks, no optimizer or OT solve."""
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import argparse,json,sys,time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    assert not args.output.exists()
    sys.path.insert(0,str(args.source.resolve()))
    import numpy as np
    import torch
    from SpaDOT.model.SpaDOT import SpaDOT
    from mechanism_core import (array_digest,fit_whitener,gp_anchor_mean,gat_anchor_mean,
        incoming_closure,all_incoming_hops,exact_influence,summarize_influence,
        coordinate_geometry_experiment,temporary_inducing_coordinates,
        _basis_jacobian_energy)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    torch.manual_seed(20261104)
    rng=np.random.default_rng(20261104)
    # Independently known linear Jacobian validates raw/whitened orientation.
    test_expression=torch.tensor(rng.normal(size=(4,3)),dtype=torch.float64,requires_grad=True)
    linear=torch.tensor(rng.normal(size=(12,10)),dtype=torch.float64)
    test_whitener,_=fit_whitener(rng.normal(size=(30,10)))
    known=_basis_jacobian_energy(test_expression.reshape(-1)@linear,test_expression,test_whitener)
    expected_raw=linear.reshape(4,3,10).square().sum((1,2)).numpy()
    expected_white=(linear@torch.tensor(test_whitener)).reshape(4,3,10).square().sum((1,2)).numpy()
    assert np.allclose(known['raw_energy'],expected_raw,atol=1e-12,rtol=1e-12)
    assert np.allclose(known['whitened_energy'],expected_white,atol=1e-12,rtol=1e-12)
    n,genes=9,7
    raw=np.column_stack([np.arange(n,dtype=float),np.arange(n,dtype=float)**2/11])
    mean,scale=raw.mean(0),raw.std(0)
    X=torch.tensor((raw-mean)/scale,dtype=torch.float64)
    Y=torch.tensor(rng.normal(size=(n,genes)),dtype=torch.float64)
    # Directed incoming chain: node0 can reach anchor5 in5hops, not3.
    edge=torch.tensor([list(range(n))+list(range(n-1)),list(range(n))+list(range(1,n))],dtype=torch.long)
    cfg=dict(input_dim=genes,z_dim=20,dtype=torch.float64,device=torch.device('cpu'),
        svgp_encoder_layers=[12,11],gat_encoder_hidden=512,gat_attention_heads=4,
        decoder_layers=[12,11],kernel_type='Gaussian',kernel_scale=.1)
    dl=dict(inducing_points={'t0':X[[0,2,4,6,8]].numpy()},N_train={'t0':n})
    model=SpaDOT(cfg,dl).to(dtype=torch.float64).eval()
    for p in model.parameters():p.requires_grad_(False)
    before={name:array_digest(value) for name,value in model.state_dict().items()}
    with torch.no_grad():
        native=model.all_latent_samples(X,Y,edge,'t0')
        direct_gp=gp_anchor_mean(model,'t0',X,Y,np.arange(n)).numpy()
        direct_gat=model.GATEncoder(Y,edge)[0].numpy()
    assert np.array_equal(np.column_stack([direct_gp,direct_gat]),native)
    anchor=5
    nodes,hops3=incoming_closure(edge,n,anchor)
    assert nodes.tolist()==[2,3,4,5] and hops3[0]==-1
    fullhops=all_incoming_hops(edge,n,anchor)
    assert fullhops[0]==5 and fullhops[8]==-1
    with torch.no_grad():
        sub,_=gat_anchor_mean(model,Y,edge,anchor,nodes)
    assert np.allclose(sub.numpy(),native[anchor,10:],atol=1e-10,rtol=1e-10)
    whiteners={};metadata={}
    for branch,z in [('GP',native[:,:10]),('GAT',native[:,10:])]:
        whiteners[branch],metadata[branch]=fit_whitener(z)
    influence=exact_influence(model,'t0',X,Y,edge,anchor,whiteners,native[anchor])
    assert influence['GAT']['raw_energy'][0]==0 and influence['GAT']['raw_energy'][8]==0
    assert influence['GP']['raw_energy'].sum()>0 and influence['GAT']['raw_energy'].sum()>0
    # Independent native fullgraph GAT derivative agrees with closure derivative.
    expression=Y.clone().requires_grad_(True)
    outputs=model.GATEncoder(expression,edge)[0][anchor]
    native_energy=torch.zeros(n,dtype=torch.float64)
    for dim in range(10):
        grad,=torch.autograd.grad(outputs[dim],expression,retain_graph=dim<9)
        native_energy+=grad.square().sum(dim=1)
    assert np.allclose(native_energy.numpy(),influence['GAT']['raw_energy'],atol=1e-11,rtol=1e-10)
    # Directional finite differences validate COMPLETE GP μ+noise derivative.
    # Test a source beyond GAT support and one in its direct support.
    checks=[]
    for branch in ['GP','GAT']:
        for source in [0,4]:
            expression=Y.clone().requires_grad_(True)
            if branch=='GP':out=gp_anchor_mean(model,'t0',X,expression,[anchor])[0]
            else:out=model.GATEncoder(expression,edge)[0][anchor]
            direction=torch.tensor(rng.normal(size=genes),dtype=torch.float64)
            direction/=direction.norm()
            vector=torch.tensor(rng.normal(size=10),dtype=torch.float64)
            vector/=vector.norm()
            gradient,=torch.autograd.grad(out@vector,expression)
            derivative=float(gradient[source]@direction)
            eps=1e-5
            delta=torch.zeros_like(Y);delta[source]=eps*direction
            with torch.no_grad():
                if branch=='GP':
                    plus=gp_anchor_mean(model,'t0',X,Y+delta,[anchor])[0]
                    minus=gp_anchor_mean(model,'t0',X,Y-delta,[anchor])[0]
                else:
                    plus=model.GATEncoder(Y+delta,edge)[0][anchor]
                    minus=model.GATEncoder(Y-delta,edge)[0][anchor]
            finite=float(((plus-minus)@vector)/(2*eps))
            assert np.isclose(derivative,finite,atol=2e-7,rtol=2e-4),(branch,source,derivative,finite)
            checks.append(dict(branch=branch,source=source,analytic=derivative,
                finite_difference=finite,absolute_error=abs(derivative-finite)))
    summary=summarize_influence(influence['GP']['whitened_energy'],X.numpy(),fullhops,anchor)
    assert np.isclose(summary['normalized_mass'].sum(),1)
    assert np.isclose(summary['normalized_nonself_mass'].sum(),1)
    assert np.isclose(sum(x['total_mass'] for x in summary['hop_rows']),1)
    null=summarize_influence(np.zeros(n),X.numpy(),fullhops,anchor)
    assert null['degenerate_zero_total'] and null['degenerate_zero_nonself']
    assert null['self_mass'] is None and null['radius90'] is None
    assert all(row['total_mass'] is None for row in null['hop_rows'])
    # Uncertainty is an active expression-dependent path, not fixed externally.
    expression=Y.clone().requires_grad_(True)
    target=gp_anchor_mean(model,'t0',X,expression,[anchor])[0].sum()
    full_gradient,=torch.autograd.grad(target,expression)
    expression=Y.clone().requires_grad_(True)
    pseudo_mean,noise=model.SVGPEncoder(expression)
    conditional=[]
    for dim in range(10):
        conditional.append(model.svgp_dict['t0'].approximate_posterior_params(
            X[[anchor]],X,pseudo_mean[:,dim],noise[:,dim].detach())[0][0])
    no_noise_gradient,=torch.autograd.grad(torch.stack(conditional).sum(),expression)
    uncertainty_gradient_difference=float((full_gradient-no_noise_gradient).abs().max())
    assert uncertainty_gradient_difference>1e-12
    # A graph callback using true directed-nearest neighbors on untied points.
    def graph_builder(points):
        distances=((points[:,None,:]-points[None,:,:])**2).sum(2)
        np.fill_diagonal(distances,np.inf)
        nearest=np.argmin(distances,axis=1)
        result=np.column_stack([np.repeat(np.arange(n),2),np.column_stack([np.arange(n),nearest]).reshape(-1)])
        result=result[np.lexsort((result[:,1],result[:,0]))]
        return torch.tensor(result.T,dtype=torch.long)
    geometry_edges=graph_builder(raw)
    with torch.no_grad():
        geometric_baseline=np.column_stack([gp_anchor_mean(model,'t0',X,Y,np.arange(n)).numpy(),
            model.GATEncoder(Y,geometry_edges)[0].numpy()])
    geometry=coordinate_geometry_experiment(model,'t0',X,Y,geometry_edges,raw,mean,scale,
        graph_builder,whiteners,geometric_baseline)
    assert all(row['rebuilt_graph_exact'] for row in geometry['rows'])
    assert geometry['inducing_restored_exactly']
    assert all(row['GAT']['max_abs_raw']==0 for row in geometry['rows'])
    identity=next(row for row in geometry['rows'] if row['transform']=='identity')
    reflection=next(row for row in geometry['rows'] if row['transform']=='rigid_reflection')
    assert identity['GP']['max_abs_raw']<1e-9 and reflection['GP']['max_abs_raw']<1e-9
    # Guaranteed rollback even when an intervention raises an exception.
    original=array_digest(model.svgp_dict['t0'].inducing_index_points)
    try:
        with temporary_inducing_coordinates(model.svgp_dict['t0'],X[:2]):
            raise RuntimeError('fixture rollback')
    except RuntimeError:pass
    assert array_digest(model.svgp_dict['t0'].inducing_index_points)==original
    assert before=={name:array_digest(value) for name,value in model.state_dict().items()}
    report=dict(all_passed=True,no_optimizer_steps=True,no_OT_solve=True,CPU_only=True,
        native_forward_exact=True,fullgraph_vs_closure_gradient_verified=True,
        complete_GP_and_GAT_finite_difference_checks=checks,whitening=metadata,
        linear_Jacobian_whitener_orientation_verified=True,
        zero_total_influence_explicitly_undefined=True,
        uncertainty_path_gradient_max_difference=uncertainty_gradient_difference,
        exact10basis_influence=True,coordinate_controls=geometry['rows'],
        inducing_rollback_verified=True,all_learned_tensors_unchanged=True,
        finished_unix=time.time())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['whitening','coordinate_controls']},indent=2))


if __name__=='__main__':main()
