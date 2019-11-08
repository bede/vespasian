import os
import warnings

from math import ceil
from pprint import pprint

import argh

from vespasian import vespasian, util, __version__



def configure_warnings(show_warnings):
    '''Show or suppress warnings, mainly for TreeSwift Tree.mrca() operations'''
    if show_warnings:
        warnings.filterwarnings('always')
    else:
        warnings.filterwarnings('ignore')


def infer_gene_trees(input: 'path to directory containing gene families',
                     tree: 'path to newick formatted species tree',
                     output: 'path to output directory' = 'gene-trees',
                     separator: 'character separating taxon name and identifier(s)' = '|',
                     warnings: 'show warnings' = False,
                     progress: 'show progress bar' = False):
    '''Create gene trees by pruning a given species tree'''
    configure_warnings(warnings)
    vespasian.infer_gene_trees(input, tree, output, separator, progress)


def codeml_setup(input: 'path to directory containing aligned gene families',
                 gene_trees: 'path to directory containing gene trees',
                 branches: 'path to yaml file containing branches to be labelled' = None,
                 output: 'path to output directory' = 'codeml',
                 separator: 'character separating taxon name and identifier(s)' = '|',
                 strict: 'label only branches with all taxa present in tree' = False,
                 threads: 'number of parallel workers' = ceil(os.cpu_count()/2),
                 warnings: 'show warnings' = False,
                 progress: 'show progress bar' = False):
    '''Create suite of branch and branch-site codeml environments'''
    configure_warnings(warnings)
    if not branches:
        print('No branch file supplied. Branch-site tests will not be configured.')
    vespasian.codeml_setup(input, gene_trees, branches, output, separator, strict, threads, progress)


def report(input: 'path to codeml_setup() output directory',
           output: 'path to output directory' = 'report-codeml'):
    '''Perform likelihood ratio tests and and report positively selected sites'''
    vespasian.report(input, output)
    print('''                            ,*//(##((##((*((,.                                 
                        *##/(###(*/(#(//(/((###%%%%#/                          
                    /(#(#((*****/**////*,,*(**/*((%%%%%%/                      
                 ///(/***,/***/**,**,,*//////**/(((###%%&%%%(.                 
              .((///,**,*//,*/*,***/****////(//((((##%%%%%%%%#%%#.             
             ((((/(///****,*,..,,,,,*,,,,,,,,,,,*////(###%%%%###%&%#           
           ###(//**,,,**,,......,.,.......,,,,,,,,,**((((#####%((###%##        
         ###(//***,,**/,,......,.,.............,,,,,,/*/((#(/((##%%###%#%.     
      ./(#(//*,,**,,,.. .....................,,,,,,,,***/(((#(/////(#((###/    
    .(((((//*,,*,,,,.......,,,,,,,.,,.........,,,,,,,****//((((*/%#((/((%(#(   
   *//((/****,,,,,,..,,,,,,,.,...,,,,,.....,,,,,,,,*,**//*/((((///*//((((/#(%  
  ,//((/*//**,,,,...,,,,,,,,,,,,,*,,...,.,,,,,***,,,,,**//((/((*(*(/##*(#/(#(/ 
 ,/(#((/***/**,,..,..,.........,......,,.,,,,,*,,,*,,***//////#*,*#*(*/&%*/%(% 
 *//#%#(*,***..,,.......,,......,.,.,,,,,,,,,,,,*****/(////(/((/////(**(///% 
.(((%%#(/*,,.,,,,..........,...,,,,,....,,,,,,,,******//((#(/(////((/#%(*(//#* 
.(##%#(*,...,,,,,.......,....,,,,..,,,,,,,,,**,,,***/(/####((/**,#%%(###*//# 
.#//%&%##(*.,,,,,,,,,,.....,,..,,,,,,,**,,,,,,,,,,,*/*(((#(((////*/###((#&@(// 
 %(&%#%#/*,,.,.,,,**,,,,,*,,,,,,,,,,,,,,*****,******(###//**//*//((*#(&%&(%% 
 .&%/&%%#%/*/*/**,,,*****,,,,,,*/**//#%%%%##/***(///#/(###(/***(//(/*#/#%%/# 
 #@/##%%(%%/,,,.,,,/&@&&(/****,,*#&@&%(//*,,,,,,,*/(#(#%#((**//(((//(&&%%#(%#/ 
 ./%/%##(##***,**,,,..#@@,...,*##%(///#&&&%##(((/(#%%%#(/(////(//##%/(*%&  
  ,#&%#/**&&&&&&(,,*#%&**/(*,,,***,,,*%%/,,,,,*(%&&%##/,////***/**/(&%,.,,*%#* 
  ,#&%#(*. %@(,..     ..***,..,*//,,,,,,*,...,,,*(///,,,*,**,,*/(//(%*.,(%&&%%/
**. (#(/,..../#,...***,.*/*,,,,/(//,,,,,/##(((##(*,,..,,,,,**/*/(##&*,#@&*(((#%
/@@#,/(*,.....,/***,. .,/(,,,,*(//**,,,.,,,,,,,,,,,,,,,,,,,*////#%%%,/#,,*(/**(
.*&@&%(*.,............,*(/,,.,*#(//**,,,,,,.,,.,,,,,*******/(/(##%#**%,.,*(*,,*
 ..*(&(**,............,//*,.,,*((#(*,.,,,,.,...,,,,,,***//((###%&%/,//**,/#/(& 
   .,/%(/*............,(*,,,,,*(%((/*,,,......,,,,,,***((###%%%&%%,**#.,*/%((/ 
    ,,%(/***,........,(#**,,**(##((((,,,,.....,,,,****/(#%%%%%&%###@@%*,*(###  
    #./%/**,,,,,.,..,,&%/***/(%%%((##,**,,,,.,,,,,,**//#%%%%%%%%##(**(#*((*#   
     %#@(*****,,,,,,,.&@@@&%@@&@@@@&/.*/**,,,,**//////(#%%%#%(##/,,*%#/*(    
     .%/****,,***,.,#@@@@@@@@&%(*...,*///****/(((/((########/(##*,.,*(#(     
     ,(#&/*(/*,,*,..,,*%&@@@@&/........,,*//***//##(/((#%#((##(#%#,*,,(%,      
      (#&(*((,,,,,..****(((/,.,,**,..,,,,,*/////(###(/(##(/#((/#&&&&&%%        
       ,//#,......,*/**,,........,*/((#*,,,*((#%%(//(#(/(%((#&@@@#*          
        .&(/%/,.....,,,**/((####(///*,,..,,,,*(##%%(/(((((##((%@@&%            
         #%(%#*,,..,**,,..,*///(/****,*,,,,,,/(#%%#(///(/#%##%&@&%.            
          %%#@%*,,****,,,.,,,,,,.,..,**,,,,,*/#%&%%(///(%%&%%&@&%%             
           @%&@&/,,,.......,..,......***,.,*/(#%%%#(/(##&&&&@@@%#*             
           ,@&&@@&(,,,..................,,,(%%&&%#(/(#%&&&&@@&(#%              
           *@@@@&@@@#*,,,,,,.,...,.,,,,***(%%&&&%###%&&@@@@@%(/#(              
           @&@@@@&@@@@@###%#%%##%%%###%&&@@@&%#%%%&@@&&&%//##&.              
          ,%##&&&&&@@@@@@&@&@&@@@@@@@@@@@@@@&&&%%%&@@&%%%/(/#%%&               
          /**/*/%&&&&&@@@@@@@@@@@@@@@@@@@@&%%#%%&@&&%#(**,/#%%&/               
          #(,*///(#&&&&&&@@&&&&&&&&&&&&&%%##%%&&%%%#(*,,**(###%,               
         //#(,,,*///#%&&&&&&&&&&&&%%%%%%%#%%%%%%%(/,,,,,,/((##&                
         *(///*,,,,////((##%%%#%%&&%%#########/*,,,,,.,,*(##%&&                
         ##/****,...,**/******,,*,,*/((///////,,,,,,,,,/(##%&&&                
        ///(/*,,,,,,,..,,,,,,,*,,,,,,,,,,,,,*,,,,,,,,**(###%&@&                
        (*,*/*,,,..,,,,....,,,*//***,,,,,,,,,,,,,,,,**/#%#%%&&&                
       .#,,,,,***,,.,,,,,,..,,,,,,,....,,,,,,,,,,,,,**#%%&&&&&%/               
        ./*,,,,**,,,,,..,,,,,,,,,,,.,,,,,,,,,,,,,***((((%%%&&%#/               
           ,.,,*,,,,.,,...,,.,,,,*,**,,,,,*,*,,,,***//(#%##%%#(,               
            ,,,,,,,,,,,...,*,,.,,..,....,,,.,,,*****((###%##((#                
               ,.,,*...,......,,,,***,,,,,,,,,*****/(######(#*                 
                .,(,,,,.,.,.,.,,,,,,,,,,,,,,,,,,,*/(####%((,                   
                   ....,,,,,.,.......,.,,,,,,,,*/(#(%###(                      
                     .,..,,,*,...,,..,**,**/***/((#%#,                         
                        .,*,.....,,,,,***/((//###/.                            
                            ,,,,,/,****///*.                                   
''')
    print(f'Report written to {output}')



###################################################################################################


def reformat_environments(input: 'path to directory containing codeml environments'):
    '''Reformat vespasian codeml environments for use with legacy codeml_reader'''
    util.reformat_environments(input)


def version():
    '''Show version'''
    print(__version__)


def main():
    argh.dispatch_commands([infer_gene_trees,
                            codeml_setup,
                            report,
                            reformat_environments,
                            version])


if __name__ == '__main__':
    main()
