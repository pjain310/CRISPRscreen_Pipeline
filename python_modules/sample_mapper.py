!#/usr/bin/python

def makeSampleMap(sample_mapf):
    '''
    Prepares a dictionary which stores information about samples and their controls.
    Args:
        sample_mapf (string) : Path to sample map file.
    Returns:
         sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
    '''
    all_lines = open(sample_mapf,'r').readlines()[1:]
    sample_map = {}
    for line in all_lines:
        try:
            condition, name, control = line.strip().split('\t')
        except:
            condition, name = line.strip().split('\t')
            control = ''
        if condition in sample_map:
            sample_map[condition]['name'].append(name)
            sample_map[condition]['control'].append(control)
        else:
            sample_map[condition] = {'name' : [name], 'control': [control]}
    return sample_map

def mageck_sampler(sample_map):
    '''
    Prepares sample mapping according to mageck-rra requirements
    Args:
         sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
    Returns:
         tab separated lists of treatments and controls for each mageck run condition
    '''

    for condition in sample_map:
        compared_samples = ",".join(sample_map[condition]['name'])
        controls = [c for c in sample_map[condition]['control'] if c]
        if len(controls)>0:
            for cl in controls[0].split(","):
                print(condition+"\n"+compared_samples+"\n"+",".join(sample_map[cl]['name']))

def mageck_mle_sampler(sample_map):
    '''
    Prepares sample mapping file according to mageck-mle requirements
    Args:
         sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
    Returns:
         design table file
    '''

    conditions=list(sample_map.keys())

    #Open design table file
    f=open("dm_mageck_mle.txt","w")

    #Write header
    f.write("Sample\tBaseline\t"+"\t".join(conditions)+"\n")

    num_conditions=len(conditions)
    i=0

    #Write condition values for each of the samples
    for condition in conditions:
        controls = [c for c in sample_map[condition]['control'] if c]
        names = sample_map[condition]['name']
        if len(controls)>0:
            for cl in controls:
                f.write(cl+"\t1"+"\t0"*num_conditions+"\n")
            for name in names:
                f.write(name+"\t1"+"\t0"*(i)+"\t1"+"\t0"*(num_conditions-i-1)+"\n")
            i+=1

    f.close()

def CB2_sampler(sample_map):
    '''
    Prepares sample mapping according to CB2 requirements
    Args:
         sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
    Returns:
         tab separated table of group and sample name
    '''

    #Create files for each condition
    for condition in sample_map:
        controls = [c for c in sample_map[condition]['control'] if c]
        names = sample_map[condition]['name']
        if len(controls)>0:
            f=open("cb2_df_"+condition+".txt","w")
            f.write("{}\t{}\n".format("group","sample_name"))
            for cl in controls[0].split(","):
                for c_name in sample_map[cl]["name"]:
                    f.write("{}\t{}\n".format("Base",c_name))
            for s_name in names:
                f.write("{}\t{}\n".format(condition,s_name))
            f.close()

def PBNPA_sampler(sample_map):
    '''
    Prepares sample mapping according to PBNPA requirements
    Args:
         sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
    Returns:
         tab separated table of group and sample name to stdout
    '''

    #Create files for each condition
    for condition in sample_map:
        controls = [c for c in sample_map[condition]['control'] if c]
        names = sample_map[condition]['name']
        print(condition,controls)
        if len(controls)>0:
            f=open("pdnpa_df_"+condition+".txt","w")
            f.write("{}\t{}\n".format("Control","Treatment"))
            for cl in controls[0].split(","):
                for c_name in sample_map[cl]["name"]:
                    for s_name in names:
                        f.write("{}\t{}\n".format(c_name,s_name))
            f.close()

CB2_sampler(makeSampleMap("sample_maps/sample_map_DLD1.txt"))
