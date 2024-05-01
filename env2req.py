import ruamel.yaml

yaml = ruamel.yaml.YAML()
data = yaml.load(open('documentation-environment.yml'))

requirements = []

for dep in data['dependencies']:
    print(dep)
    if isinstance(dep, str):
        package_info = dep.split('=')
        package = package_info[0]
        try:
            package_version = package_info[1]
        except:
            package_version = None
        
        try:
            python_version = package_info[2]package = package_info-
        except: 
            python_version = None
        
        !!!!!!!!!!!!!!!!!!!!!
        requirements.append(package + '==' + package_version)
    elif isinstance(dep, dict):
        for preq in dep.get('pip', []):
            requirements.append(preq)

with open('requirements.txt', 'w') as fp:
    for requirement in requirements:
        print(requirement, file=fp)
