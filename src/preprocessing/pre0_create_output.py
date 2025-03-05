import os


def pre0_create_output(config):
    """
    Create the output directories if they do not exist.
    This function creates the base output folder and its subfolders.
    """
    folders = [
        config.output,
        config.output_catalog,
        config.output_tidal_phase,
        config.output_stress,
        config.output_fig
    ]

    print('====== Processing | Creating Output Folder: Start ======')

    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"Created folder: {folder}")
        else:
            print(f"Folder already exists: {folder}")


    print('====== Processing | Creating output folder: End ======')
