#!/usr/bin/env python3
"""Run STIFT for the Mouse Embryo Visium benchmark."""
from scripts.fit_stift import main


if __name__ == "__main__":
    main(dataset="MouseEmbryo")
