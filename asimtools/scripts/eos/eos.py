'''.Xauthority'''

from typing import Dict, Tuple
from asimtools.job import ChainedJob

def eos(
    image: dict,
    singlepoint_env_id: str,
    calc_id: str,
    preprocess_env_id: str,
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
    postprocess: bool = True,
) -> Dict:
    ''' EOS calculation '''

    sim_input = {
        'step-0': {
            'script': 'eos.preprocess',
            'env_id': preprocess_env_id,
            'args': {
                'image': image,
                'nimages': nimages,
                'scale_range': scale_range,
            }
        },
        'step-1': {
            'script': 'image_array',
            'env_id': preprocess_env_id,
            'args': {
                'images': {'image_file': 'step-0/preprocess_images_output.xyz'},
                'subscript_input': {
                    'script': 'singlepoint',
                    'env_id': singlepoint_env_id,
                    'args': {
                        'calc_id': calc_id,
                        'properties': ['energy'],
                    }
                }
            }
        },
        'step-2': {
            'script': 'eos.postprocess',
            'env_id': preprocess_env_id,
            'submit': postprocess, #Fails if previous step is in a slurm queue
            'args': {
                'images': {'pattern': '../step-1/id-*/image_output.xyz'},
            }
        }
    }

    chain = ChainedJob(sim_input, env_input=None, calc_input=None)
    chain.submit()
    return chain.get_last_output()
