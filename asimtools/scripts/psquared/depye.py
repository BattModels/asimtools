'''.Xauthority'''

from typing import Dict, Tuple
from asimtools.job import ChainedJob, branch

@branch
def eos(
    config_input: Dict,
    image: dict,
    preprocess_config_id: str,
    singlepoint_config_id: str,
    nimages: int = 5,
    scale_range: Tuple[float] = (0.95, 1.05),
    **kwargs
) -> Tuple[list,Dict]:
    ''' EOS calculation '''

    sim_input = {
        'step-0': {
            'script': 'eos.singlepoint_calculations',
            'config_id': preprocess_config_id,
            'args': {
                'image': image,
                'singlepoint_config_id': singlepoint_config_id,
                'preprocess_config_id': preprocess_config_id,
                'nimages': nimages,
                'scale_range': scale_range,
            },
        },
        'step-1': {
            'script': 'eos.postprocess',
            'config_id': preprocess_config_id,
            'args': {
                'step0_dir': '../step-0',
                'scale_range': scale_range,
            },
        },
        'step-2': {
            'script': 'depye.gthermo_calculations',
            'config_id': preprocess_config_id,
            'args': {
                'step1_dir': '../step-1',
                'scale_range': scale_range,
            },
        },
    }
    chain = ChainedJob(config_input, sim_input)
    job_ids = chain.submit()
    return job_ids, chain.get_last_output()
