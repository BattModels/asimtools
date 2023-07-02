from pathlib import Path


src_dir = Path("./")
for file in src_dir.iterdir():
    if str(file).endswith('.rst'):
        print("Processed RST file:", file)
        with open(file, "r") as f:
            lines = f.read()

        junk_strs = ["Submodules\n----------", "Subpackages\n-----------"]

        for junk in junk_strs:
            lines = lines.replace(junk, "")

        lines = lines.replace(" module\n=", "\n")

        with open(file, "w") as f:
            f.write(lines)
