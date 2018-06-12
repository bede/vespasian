import os
import sys

family_paths = [fn for fn in os.listdir(root_path) if not fn.startswith('.')]
families_paths = {fn.partition('|')[0]: fn for fn in families_paths}

models_presets = {
    'm0': {'model': 0, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 1},
    '2ratios': {'model': 2, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 2},
    'm1Neutral': {'model': 0, 'NSsites': 1, 'fix_omega': 0, 'ncatG': 2},
    'm2Selection': {'model': 0, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 3},
    'm3Discrtk2': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 2},
    'm3Discrtk3': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 3},
    'm7': {'model': 0, 'NSsites': 7, 'fix_omega': 0, 'ncatG': 10},
    'm8': {'model': 0, 'NSsites': 8, 'fix_omega': 0, 'ncatG': 10},
    'm8a': {'model': 0, 'NSsites': 8, 'fix_omega': 1, 'ncatG': 10, 'omega': 1},
    'modelA': {'model': 2, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 4},
    'modelAnull': {'model': 2, 'NSsites': 2, 'fix_omega': 1, 'ncatG': 4, 'omega': 1},
}

class ControlFile():
    '''Create a CodeML control (.ctl) file'''
    def __init__(self, model_name, model=0, NSsites=8, fix_omega=1, ncatG=10, omega=1):
        self.seqfile = 'align.phy'
        self.treefile = 'tree'
        self.outfile = 'out'
        self.noisy = 3
        self.verbose = 0
        self.runmode = 0
        self.seqtype = 1
        self.CodonFreq = 2
        self.aaDist = 0
        self.aaRatefile = 'wag.dat'
        self.model = model
        self.NSsites = NSsites
        self.icode = 0
        self.fix_kappa = 0
        self.kappa = 3
        self.fix_omega = fix_omega
        self.omega = omega
        self.fix_alpha = 1
        self.alpha = 0
        self.Malpha = 0
        self.ncatG = ncatG
        self.clock = 0
        self.getSE = 0
        self.RateAncestor = 0
        self.Small_Diff = .5e-6
    
    def __str__(self):
        return self.__dict__()

    def write(self, path):
        '''Write valid .ctl file to specified path'''
        with open(path, 'w') as ctl_fh:
            pass 




def workspace_setup(family_path):
    '''Create a CodeML workspace for performing standard tests on an orthologous family'''
    
    ctl_defaults = {
        seqfile: 'align.phy', 
        treefile: 'tree', 
        outfile: 'out', 
        noisy: 3,
        verbose: 0,
        runmode: 0,
        seqtype: 1,
        CodonFreq: 2,
        aaDist: 0,
        aaRatefile: 'wag.dat',
        model: 0,
        NSsites: 0,
        icode: 0,
        fix_kappa: 0,
        kappa: 3,
        fix_omega: 0,
        omega: 1,
        fix_alpha: 1,
        alpha: 0,
        Malpha: 0,
        ncatG: 1,
        clock: 0,
        getSE: 0,
        RateAncestor: 0,
        Small_Diff: .5e-6
    }

    models_presets = {
        'm0': {'model': 0, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 1},
        '2ratios': {'model': 2, 'NSsites': 0, 'fix_omega': 0, 'ncatG': 2},
        'm1Neutral': {'model': 0, 'NSsites': 1, 'fix_omega': 0, 'ncatG': 2},
        'm2Selection': {'model': 0, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 3},
        'm3Discrtk2': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 2},
        'm3Discrtk3': {'model': 0, 'NSsites': 3, 'fix_omega': 0, 'ncatG': 3},
        'm7': {'model': 0, 'NSsites': 7, 'fix_omega': 0, 'ncatG': 10},
        'm8': {'model': 0, 'NSsites': 8, 'fix_omega': 0, 'ncatG': 10},
        'm8a': {'model': 0, 'NSsites': 8, 'fix_omega': 1, 'ncatG': 10, 'omega': 1},
        'modelA': {'model': 2, 'NSsites': 2, 'fix_omega': 0, 'ncatG': 4},
        'modelAnull': {'model': 2, 'NSsites': 2, 'fix_omega': 1, 'ncatG': 4, 'omega': 1},
    }

    models_with_fixed_omega = set('m8a', 'modelAnull')  # Models with fixed omega value of 1

    omega_values = [0, 1, 2, 10]
    
    model_dirs = [f'{family_path}/{m}' for m in models_presets.keys()]

    for model, presets in models_presets.items():
        os.mak
        ctl_config = ctl_defaults**
        ctl_config.update(presets)
        if model not in models_with_fixed_omega:
            for omega in omega_values:
                    ctl_config['omega'] = omega


for model, presets in models_presets.items():
    if model not in models_with_fixed_omega:
        for omega in omega_values:
            ctl = ControFile(model, presets*)



















class CodemlWorkspace:
    def __init__(self, preset):
        self.id = 


    def write():
        pass



class ControlFile():
    def __init__(self, preset):
        self.preset = preset

    def write():
        pass
