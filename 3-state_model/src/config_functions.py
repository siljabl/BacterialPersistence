def save_config(config, dir):
    '''
    Saving simulation parameters
    config: dict of parameters
    dir: where config is saved
    '''
    with open(f"{dir}/config.txt", 'w') as f: 
        #f.write('Units\nlenght:pixels\ntime:frames\n\n') 

        for key, value in config.items():  
            f.write('%s:%s\n' % (key, value))


def read_config(dir):
    '''
    Importing simulation parameters
    config: dict of parameters
    dir: where config is saved
    '''   
    with open(f"{dir}/config.txt", 'r') as file:
        lines = file.readlines()

    config = {}    
    for line in lines:
        key, value = line.strip().split(':')
        config[key.strip()] = value.strip()

    return config