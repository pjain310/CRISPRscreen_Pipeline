#!/usr/bin/python

from abc import ABC,abstractmethod


#Defining base class which will be inherited by all tool classes - contains base function for mapping samplemap file to dict
class sampler(ABC):
    '''Base class for all software runs with makeSampleMapDict and SampleMapper functions.'''

    def __init__(self, sample_mapf, sample_map, *args, **kwargs):
        super().__init__()
        self.sample_mapf= sample_mapf
        self.sample_map=sample_map

    def makeSampleMapDict(self):
        '''
        Prepares a dictionary which stores information about samples and their controls.
        Args:
            sample_mapf (string) : Path to sample map file.
        Returns:
             sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
        '''
        all_lines = open(self.sample_mapf,'r').readlines()[1:]
        self.sample_map = {}
        for line in all_lines:
            try:
                condition, name, control = line.strip().split('\t')
            except:
                condition, name = line.strip().split('\t')
                control = ''
            if condition in self.sample_map:
                self.sample_map[condition]['name'].append(name)
                self.sample_map[condition]['control'].append(control)
            else:
                self.sample_map[condition] = {'name' : [name], 'control': [control]}


        return self.sample_map

    def SampleMapper(self):
        '''Placeholder function to be inherited by each class'''
        raise NotImplementedError()

#Defining class for mageckrra - contains SampleMapper specific to tool
class mageckrra(sampler):

    def __init__(self, sample_mapf, sample_map, *args, **kwargs):
        super().__init__(sample_mapf, sample_map)


    def SampleMapper(self):
        '''
        Prepares sample mapping according to mageck-rra requirements
        Args:
             sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
        Returns:
             print to stdout: tab separated lists of treatments and controls for each mageck run condition
        '''

        for condition in self.sample_map:
            compared_samples = ",".join(self.sample_map[condition]['name'])
            controls = [c for c in self.sample_map[condition]['control'] if c]
            if len(controls)>0:
                for cl in controls[0].split(","):
                    print(condition+"\n"+cl+"\n"+compared_samples+"\n"+",".join(self.sample_map[cl]['name']))

#Defining class for mageckmle - contains SampleMapper specific to tool
class mageckmle(sampler):

    def __init__(self, sample_mapf, sample_map, *args, **kwargs):
        super().__init__(sample_mapf, sample_map)

    def SampleMapper(self):
        '''
        Prepares sample mapping file according to mageck-mle requirements
        Args:
             sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
        Returns:
             design table file
        '''
        conditions=list(self.sample_map.keys())

        #Open design table file
        f=open("temp/sample_maps/mle_df__conditions_vs_control.txt","w")

        #Write header
        f.write("Sample\t"+"\t".join(conditions)+"\n")

        num_conditions=len(conditions)-1 #Assuming only one control per condition (naive)
        i=0

        #Write condition values for each of the samples
        for condition in conditions:
            controls = [c for c in self.sample_map[condition]['control'] if c]
            names = self.sample_map[condition]['name']
            if len(controls)>0:
                for cl in self.sample_map[controls[0]]['name']:
                    f.write(cl+"\t1"+"\t0"*num_conditions+"\n")
                for name in names:
                    f.write(name+"\t1"+"\t0"*(i)+"\t1"+"\t0"*(num_conditions-i-1)+"\n")
                i+=1

        f.close()

        print("temp/sample_maps/mle_df__conditions_vs_control.txt")
        print(controls[0])

#Defining class for cb2 - contains SampleMapper specific to tool
class cb2(sampler):

    def __init__(self, sample_mapf, sample_map, *args, **kwargs):
        super().__init__(sample_mapf, sample_map)

    def SampleMapper(self):
        '''
        Prepares sample mapping according to CB2 requirements
        Args:
             sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
        Returns:
             tab separated table of group and sample name
        '''

        #Create files for each condition
        for condition in self.sample_map:
            controls = [c for c in self.sample_map[condition]['control'] if c]
            names = self.sample_map[condition]['name']
            if len(controls)>0:
                for cl in controls[0].split(","):
                    f=open("temp/sample_maps/cb2_df__"+condition+"_vs_"+cl+".txt","w")
                    f.write("{}\t{}\n".format("group","sample_name"))
                    for c_name in self.sample_map[cl]["name"]:
                        f.write("{}\t{}\n".format(cl,c_name))
                    for s_name in names:
                        f.write("{}\t{}\n".format(condition,s_name))
                    f.close()
                    print("temp/sample_maps/cb2_df__"+condition+"_vs_"+cl+".txt")


#Defining class for pbnpa - contains SampleMapper specific to tool
class pbnpa(sampler):

    def __init__(self, sample_mapf, sample_map, *args, **kwargs):
        super().__init__(sample_mapf, sample_map)

    def SampleMapper(self):
        '''
        Prepares sample mapping according to PBNPA requirements
        Args:
             sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
        Returns:
             tab separated table of group and sample name to stdout
        '''

        #Create files for each condition
        for condition in self.sample_map:
            controls = [c for c in self.sample_map[condition]['control'] if c]
            names = self.sample_map[condition]['name']
            if len(controls)>0:
                for cl in controls[0].split(","):
                    f=open("temp/sample_maps/pbnpa_df__"+condition+"_vs_"+cl+".txt","w")
                    f.write("{}\t{}\n".format("Control","Treatment"))
                    for c_name in self.sample_map[cl]["name"]:
                        for s_name in names:
                            f.write("{}\t{}\n".format(c_name,s_name))
                    f.close()
                    print("temp/sample_maps/pbnpa_df__"+condition+"_vs_"+cl+".txt")

#Defining class for bagel - contains SampleMapper specific to tool
class bagel(sampler):

    #Adding counts file to list of required arguments for this class (to get column indices)
    def __init__(self, sample_mapf, sample_map, counts_file,*args, **kwargs):
        super().__init__(sample_mapf, sample_map)
        self.counts_file=counts_file

    def SampleMapper(self):
        '''
        Prepares sample mapping according to bagel requirements
        Args:
             sample_map (dictionary) : Dictionary whose key is the condition and value is a dictionary containing the list of samples belonging to the condition and control group.
        Returns:
             print to stdout: tab separated lists of treatments and controls for each bagel run condition
        '''

        header = open(self.counts_file,"r").readline().strip().split("\t")

        for condition in self.sample_map:
            treatment_indices = []
            compared_samples = self.sample_map[condition]['name']
            controls = [c for c in self.sample_map[condition]['control'] if c]
            if len(controls)>0:
                for cl in controls[0].split(","):
                    control_indices = []
                    counts_header=header.copy()
                    for cl_sample in self.sample_map[cl]['name']:
                        control_indices.append(str(header.index(cl_sample)-1))
                        counts_header.remove(cl_sample)
                    for i in compared_samples:
                        treatment_indices.append(str(counts_header.index(i)-1))
                    cl_index=",".join(control_indices)
                    trt_index = ",".join(treatment_indices)
                    print("{}\n{}\n{}".format(condition+"_vs_"+cl,cl_index,trt_index))

if __name__ == '__main__':
    import sys
    tool = sys.argv[1]
    if tool == 'mageckrra':
        sampler = mageckrra(*sys.argv[2:],'')
    elif tool == 'mageckmle':
        sampler = mageckmle(*sys.argv[2:],'')
    elif tool == 'pbnpa':
        sampler = pbnpa(*sys.argv[2:],'')
    elif tool == 'cb2':
        sampler = cb2(*sys.argv[2:],'')
    sampler.makeSampleMapDict()
    sampler.SampleMapper()

# tester=pbnpa("../inputs/sample_maps/sample_map_DLD1.txt","")
# tester.makeSampleMapDict()
# tester.SampleMapper()
