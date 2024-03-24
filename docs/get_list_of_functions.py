#%%
import os

import yaml
import inspect
import importlib.util

def get_functions_from_file(file_path):
    spec = importlib.util.spec_from_file_location("module.name", file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    functions = inspect.getmembers(module, inspect.isfunction)
    return [func[0] for func in functions]
def update_quarto_yml(directory, python_files, quarto_yml_path):
    # Load the existing _quarto.yml file
    with open(quarto_yml_path, 'r') as file:
        quarto_yml = yaml.safe_load(file)

    # For each Python file, create a new section in the _quarto.yml file
    for python_file in python_files:
        module_name = python_file[:-3]  # Remove the .py extension
        functions = get_functions_from_file(os.path.join(directory, python_file))

        section = {
            'title': module_name.capitalize(),
            'package': f'mofography.{module_name}',
            'desc': 'Description of this module.',
            'contents': functions
        }

        quarto_yml['quartodoc']['sections'].append(section)

    # Write the updated _quarto.yml data to a new file named quarto_temp.yml
    with open('quarto_temp.yml', 'w') as file:
        yaml.safe_dump(quarto_yml, file, sort_keys=False)

    return quarto_yml

# Use the function
directory = "../mofography"
quarto_yml_path = "./_quarto.yml"
python_files = ['compute.py', 'visualize.py', 'analyze.py', 'data.py', 'helper.py']
quarto_yml = update_quarto_yml(directory, python_files, quarto_yml_path)

# %%
