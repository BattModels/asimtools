'''
Collects images into a single file

author: mkphuthi@github.com
'''

from typing import Dict, Optional, Sequence
import numpy as np
from ase.io import write
from asimtools.utils import (
    get_atoms, get_images, new_db
)


def collect_images(
    images: Dict,
    out_format: str = 'xyz',
    prefixes: Sequence[str] = ['output_images'],
    splits: Optional[Sequence[float]] = (1,),
    shuffle: bool = True,
) -> Dict:
    """Collects images into one file/database and can split them into 
    multiple files/databases for ML tasks

    :param images: Images specification, see :func:`asimtools.utils.get_images`
    :type images: Dict
    :param out_format: output file format, defaults to 'extxyz'
    :type out_format: str, optional
    :param prefix: file name without extension, defaults to 'output_images'
    :type prefix: str, optional
    :param splits: Ratios to split data into, defaults to None
    :type splits: Optional[Sequence[float]], optional
    :param shuffle: shuffle images before splitting, defaults to True
    :type shuffle: bool, optional
    :return: results
    :rtype: Dict
    
    """
    if prefixes == (1):
        prefixes = [f'{prefix}-{i:03d}' for i in range(len(splits))]
    assert len(splits) == len(prefixes), \
        'Number of splits must match number of prefixes'
    images = get_images(**images)
    if shuffle:
        np.random.shuffle(images)
    datasets = []
    start_index = 0
    index_ranges = []
    for i, split in enumerate(splits):
        index_range = int(split / np.sum(splits) * len(images))
        index_ranges.append(index_range)
        dataset = images[start_index:start_index + index_range]
        if format == 'db':
            with new_db(f'{prefixes[i]}.db') as db:
                db.write(dataset)
        else:
            write(f'{prefixes[i]}.{out_format}', dataset, format=out_format)
        start_index += index_range

    return {"split": index_ranges}
