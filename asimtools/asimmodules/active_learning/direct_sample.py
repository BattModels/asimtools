from pathlib import Path
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from ase.io import write
import numpy as np
from asimtools.utils import (
    get_images,
)
from pymatgen.io.ase import AseAtomsAdaptor
from maml.sampling.direct import (
    BirchClustering, DIRECTSampler, SelectKFromClusters
)

def plot_PCAfeature_coverage(all_features, selected_indexes, method="DIRECT"):
    fig, ax = plt.subplots(figsize=(5, 5))
    selected_features = all_features[selected_indexes]
    ax.plot(
        all_features[:, 0],
        all_features[:, 1],
        "*",
        alpha=0.5,
        label=f"All {len(all_features):,} structures"
    )
    ax.plot(
        selected_features[:, 0],
        selected_features[:, 1],
        "*",
        alpha=0.5,
        label=f"{method} sampled {len(selected_features):,}",
    )
    legend = ax.legend(
        frameon=False,
        fontsize=14,
        loc="upper left",
        bbox_to_anchor=(-0.02, 1.02),
        reverse=True
    )
    for lh in legend.legendHandles:
        lh.set_alpha(1)
    ax.set_ylabel("PC 2", size=20)
    ax.set_xlabel("PC 1", size=20)

def direct_sample(
    images: dict,
    n: int | None = None,
    threshold_init: float = 0.2,
    k: int = 1,
) -> dict:

    images = get_images(**images)
    structs = [AseAtomsAdaptor.get_structure(image) for image in images]
    num_images = len(images)
    print(f'Num_images: {num_images}')
    DIRECT_sampler = DIRECTSampler(
        clustering=BirchClustering(n=n, threshold_init=threshold_init),
        select_k_from_clusters=SelectKFromClusters(k=k),
        structure_encoder='MatGLM3GNet',
        # model_path='/Users/keith/repos/matpes/weights_PBE'
    )
    print('Fitting...')
    selection = DIRECT_sampler.fit_transform(structs)

    stats = {
        'num_images': len(images),
        'num_selected': len(selection),
    }

    DIRECT_selection = selection
    explained_variance = DIRECT_sampler.pca.pca.explained_variance_ratio_

    npca = 10
    DIRECT_selection["PCAfeatures_unweighted"] = \
        DIRECT_selection["PCAfeatures"] / explained_variance[:npca]

    fig, ax = plt.subplots()
    ax.plot(
        range(1, npca+1),
        explained_variance[:npca] * 100,
        "o-",
    )
    ax.set_xlabel(r"i$^{\mathrm{th}}$ PC", size=20)
    ax.set_ylabel("Explained variance", size=20)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    all_features = DIRECT_selection["PCAfeatures_unweighted"]
    selected_indexes = DIRECT_selection["selected_indexes"][:]
    plot_PCAfeature_coverage(all_features, selected_indexes)

    selected_images = [images[i] for i in selected_indexes]
    write('selected.xyz', selected_images)
    write('all.xyz', images)

    with open('selected_indices.dat', 'w') as fp:
        fp.write(','.join([str(s) for s in selected_indexes]))

    results = {
        'num_images': len(images),
        'num_selected': len(selection),
    }

    return results