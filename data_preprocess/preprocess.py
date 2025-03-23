import sys
import scanpy as sc
from utils import LangCellTranscriptomeTokenizer

if len(sys.argv) != 3:
    print("Usage: python script.py <input_adata.h5ad> <output_tokenized_dataset>")
    sys.exit(1)

input_adata_path = sys.argv[1]
output_dataset_path = sys.argv[2]


data = sc.read_h5ad(input_adata_path)
data.obs['n_counts'] = data.X.sum(axis=1)
print("here1")
data.var['ensembl_id'] = data.var.index
print("here2")

tk = LangCellTranscriptomeTokenizer(dict([(k, k) for k in data.obs.keys()]), nproc=4)
tokenized_cells, cell_metadata = tk.tokenize_anndata(data)
tokenized_dataset = tk.create_dataset(tokenized_cells, cell_metadata)

tokenized_dataset.save_to_disk(output_dataset_path)

print(f"Tokenized dataset saved to: {output_dataset_path}")


