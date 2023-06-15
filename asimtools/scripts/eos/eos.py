'''.Xauthority'''

from typing import Dict, Tuple
from asimtools.job import ChainedJob, branch

# @branch
def eos(
    image: dict,
    singlepoint_env_id: str,
    singlepoint_calc_id: str,
    preprocess_env_id: str,
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
) -> Tuple[list,Dict]:
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
                        'calc_id': singlepoint_calc_id,
                        'properties': ['energy'],
                    }
                }
            }
        },
        'step-2': {
            'script': 'eos.postprocess',
            'env_id': preprocess_env_id,
            'args': {
                'images': {'pattern': '../step-1/id-*/image_output.xyz'},
                'scale_range': scale_range,
            }
        }
    }
    # preprocess_step = UnitJob(sim_input={
    #     'script': 'eos.preprocess',
    #     'env_id': preprocess_env_id,
    #     'args': {
    #         'image': image,
    #         'nimages': nimages,
    #         'scale_range': scale_range,
    # }})

    

    chain = ChainedJob(sim_input, env_input=None, calc_input=None)
    job_ids = chain.submit()
    return chain.get_last_output()
