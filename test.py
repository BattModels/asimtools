from pathlib import Path
version_dict = {}
with open(Path(__file__).parents[0] / 'asimtools/_version.py') as fp:
    exec(fp.read(), version_dict)

print(version_dict['__version__'])
