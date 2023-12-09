#!/usr/bin/env python
'''
asimmodule for unit tests and debugging, does nothing for specified duration
'''
from time import sleep
from typing import Dict, Optional

def do_nothing(duration: Optional[int] = 5) -> Dict:
    """Sleep for the specified duration

    :param duration: time in seconds, defaults to 60
    :type duration: int, optional
    :return: Dictionary with duration slept for
    :rtype: Dict
    """
    sleep(duration)
    results = {'duration': duration}
    return results
