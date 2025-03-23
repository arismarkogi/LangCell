import sys
import scanpy as sc
from utils import LangCellTranscriptomeTokenizer

if len(sys.argv) != 3:
    print("Usage: python script.py <input_adata.h5ad> <output_tokenized_dataset>")
    sys.exit(1)

input_adata_path = sys.argv[1]
output_dataset_path = sys.argv[2]

try:
    data = sc.read_h5ad(input_adata_path)
    data.obs['n_counts'] = data.X.sum(axis=1)
    data.var['ensembl_id'] = data.var['feature_id']

    tk = LangCellTranscriptomeTokenizer(dict([(k, k) for k in data.obs.keys()]), nproc=4)
    tokenized_cells, cell_metadata = tk.tokenize_anndata(data)
    tokenized_dataset = tk.create_dataset(tokenized_cells, cell_metadata)

    tokenized_dataset.save_to_disk(output_dataset_path)

    print(f"Tokenized dataset saved to: {output_dataset_path}")

except FileNotFoundError:
    print(f"Error: Input file not found: {input_adata_path}")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
