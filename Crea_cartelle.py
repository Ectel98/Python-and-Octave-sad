from pathlib import Path

for i in range (15,21):
    Path('intervalli risultanti a' + str(i)).mkdir(parents=True, exist_ok=True)
