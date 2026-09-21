import ast

with open('/home/jkniblo/fusions/multidomain_metapredictv3_Dec17.csv', 'r') as f:
    for i, line in enumerate(f):
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        try:
            ast.literal_eval(parts[1])
        except Exception as e:
            print(f"Line {i+1}: {line}")
            print(f"  Error: {e}")
