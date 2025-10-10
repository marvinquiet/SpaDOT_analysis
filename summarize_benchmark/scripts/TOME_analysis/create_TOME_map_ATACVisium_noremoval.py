import json
import os, sys

time_point = ["t0", "t1", "t2"]
time_point_n = len(time_point)
time_point_id = {}
for i in range(0,time_point_n):
    time_point_id[time_point[i]] = i


### read edge info 
edge = {}; node_all = []; node_each = {}; main_edge = set()
base_dir = "/media/sda/wenjinma"
result_dir = base_dir+os.sep+"projects/spatialATAC-explore"
downstream_dir = result_dir+os.sep+"summarize_benchmark/ATAC-Visium_noremoval_figures/ATAC-Visium_noremoval_downstream"+os.sep+'TOME'
file = open(downstream_dir+os.sep+"/mouse_edge.txt")
for line in file:
    l = line.rstrip().split('\t')
    if l[0] not in node_all:
        node_all.append(l[0])
    if l[1] not in node_all:
        node_all.append(l[1])
    edge[l[0]] = edge.get(l[0], [])
    edge[l[0]].append(l[1])
    num = time_point_id[l[0].split(':')[0]]
    node_each[num] = node_each.get(num, [])
    if l[0] not in node_each[num]:
        node_each[num].append(l[0])
    main_edge.add((l[0],l[1]))

### read edge weight file
main_weight = {}; extra_weight = {}
file = open(downstream_dir+os.sep+"/edge_prob.txt")
for line in file: 
    l = line.rstrip().split('\t')
    if float(l[2])>0.7:
        xx = "specific"
    elif float(l[2])>0.3:
        xx = "additional"
    if (l[0],l[1]) in main_edge:
        if float(l[2]) > 0.3:
            main_weight[l[1]] = float(l[2])
        else:
            main_weight[l[1]] = 0
    else:
        if float(l[2]) > 0.3:
            extra_weight[(l[0],l[1])] = float(l[2])


### read group information, used for the color of the node
file = open(downstream_dir+os.sep+"/celltype_group.txt")
i = 1
coor = {}
node_group = {}
for line in file:
    l = line.rstrip().split(" ")
    node_group[l[0]] = int(l[1])
    coor[l[0]] = i
    i += 1
file.close()

### create the info for all the node
dat = {}
for i in node_all:
    dat[i] = {'name':i}
    if i in edge:
        dat[i]['children'] = []
    if i in main_weight:
        dat[i]['edge_weight'] = main_weight[i]
    else:
        print(i)
    ### 10. coors of x (used to create the plot, the cell type order)
    if i in coor:
        dat[i]['fx'] = str(coor[i])
    else:
        print(i)
    ### 11. group of celltypes. 1: neuroectoderm; 2: surface ectoderm; 3: endoderm; 4: mesoderm; 5: exe germ layer; 0: other
    ### 1: green, 2: blue, 3: red; 4: yellow; 5: brown; 0: black
    color_map = {
        1: "#F8766D", 
        2: "#B79F00",
        3: "#00BA38",
        4: "#4EAEF0",
        5: "#D772EC"}
        # "t0:1":"#F8766D", "t0:4":"#B79F00", "t0:3":"#00BA38", "t0:2":"#4EAEF0", "t0:0":"#D772EC", "t1:0":"#4EAEF0", "t1:1":"#F8766D", "t1:2":"#00BA38", "t1:3":"#D772EC", "t1:4":"#B79F00", "t2:1":"#F8766D", "t2:3":"#4EAEF0", "t2:4":"#B79F00", "t2:0":"#990099", "t2:2":"#FF99FF"}
    if i in node_group:
        dat[i]['node_group'] = color_map[node_group[i]]
    else:
        print(i)


### add extra node first
for i in extra_weight:
    tmp = {'name': i[1], 'edge_weight': extra_weight[i], 'fx': str(coor[i[1]]), 'node_group': color_map[node_group[i[1]]]}
    if i[0] in dat:
        dat[i[0]]['children'] = dat[i[0]].get('children',[])
        dat[i[0]]['children'].append(tmp)
        print(dat[i[0]]['children'])


### connect each time point
# for i in range(time_point_n-2, 0, -1):
for i in range(time_point_n-1):
    for j in node_each[i]:
        for k in edge[j]:
            dat[j]['children'].append(dat[k])

dat_json = {'name': "node", "children": [], "node_group": "#000000"}
for key, value in dat.items():
    timepoint = key.split(':')[0]
    if timepoint == 't0':
        dat_json['children'].append(value)

with open(downstream_dir+os.sep+"ATACVisium_noremoval.json", 'w') as json_file:
    json.dump(dat_json, json_file)