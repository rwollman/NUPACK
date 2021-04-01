
import os, json
from tempfile import TemporaryDirectory
from pathlib import Path

def run_notebook(path):
    import nbformat
    from nbconvert.preprocessors import ExecutePreprocessor

    with path.open('r') as f:
        nb = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    ep.preprocess(nb)



def test_notebooks():
    env = os.environ.get('NUPACK_NOTEBOOK_TEST_DIR')
    if env is None:
        return
    env = Path(env)

    for notebook in env.glob('**/*.ipynb'):
        if '.ipynb_checkpoints' in str(notebook):
            continue
        try:
            run_notebook(notebook)
        except Exception as e:
            raise RuntimeError('Notebook %s failed' % notebook) from e

