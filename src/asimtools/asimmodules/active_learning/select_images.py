from typing import Dict, Sequence, Optional
import os
import random
import pandas as pd
from ase.io import write
from asimtools.utils import (
    get_images,
)

def select_images(
    image_file: os.PathLike,
    deviation_csv: os.PathLike,
    thresholds: Sequence[float],
    column: Optional[str] = 'force_mean_std',
    max_images: Optional[int] = None,
) -> Dict:
    """Select images based on deviation values

    :param image_file: Path to image_file with images to select from
    :type image_file: os.PathLike
    :param deviation_csv: Path to csv file from compute_deviation, row indices
        should match indices in image_file
    :type deviation_csv: os.PathLike
    :param thresholds: Thresholds for selection provided as (min,max)
    :type thresholds: Sequence[float]
    :param column: Column of csv to use, defaults to 'force_mean_std'
    :type column: str, optional
    :param max_images: Maximum number of images to select, defaults to None
    :type max_images: int, optional
    :return: Dictionary with stats on selected images
    :rtype: Dict
    """
    
    dev_df = pd.read_csv(deviation_csv, sep=',')
    selected_idxs = dev_df.index[
        (dev_df[column] > thresholds[0]) & (dev_df[column] < thresholds[1])
    ]

    images = get_images(image_file=image_file)
    selected_images = [images[i] for i in selected_idxs]
    random.shuffle(selected_images)
    num_selected = len(selected_images)
    if max_images is not None:
        selected_images = selected_images[:max_images]
    write('selected_images.xyz', selected_images)
    results = {
        'num_selected': num_selected,
        'fraction_selected': num_selected/len(images),
        'final_size': len(selected_images),
        'selected_idxs': list(selected_idxs),
    }
    return results
    