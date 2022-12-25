def make_raspa_input_for_cif(folder_path, raspa_settings_dict, other_files_to_copy, templates_dir, template_name):
    import jinja2 as jj 
    import shutil
    import os 
    ej = jj.Environment(loader = jj.FileSystemLoader(templates_dir))
    tlt = ej.get_template(template_name)
    simulation_input_string = tlt.render(**raspa_settings_dict)
    with open(os.path.join(os.path.abspath(folder_path), 'simulation.input'), 'w') as f: 
        f.write(simulation_input_string)
    for f in other_files_to_copy:
        shutil.copy(os.path.abspath(f), os.path.abspath(folder_path))