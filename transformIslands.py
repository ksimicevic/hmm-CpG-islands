import pandas as pd

path = input("File path:")
res_folder = input("Folder to save transformed file: ")

df = pd.read_csv(path)
df[['chromStart', 'chromEnd']].to_csv(res_folder + "result.csv", index=False, header=False)
print("Done")