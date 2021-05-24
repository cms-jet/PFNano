import yaml
import argparse
import os
import pprint
import re
import copy

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor coffea files')
parser.add_argument('-c', '--card', '--yaml', dest='card', default='card_example.yml', help='Crab yaml card')
parser.add_argument('--make', action='store_true', help="Make crab configs according to the spec.")
parser.add_argument('--submit', action='store_true', help="Submit configs created by ``--make``.")
parser.add_argument('--status', action='store_true', help="")
parser.add_argument("--test", type=str2bool, default='True', choices={True, False}, help="Test submit - only 1 file, don't publish.")
args = parser.parse_args()

with open(args.card, 'r') as f:
    card = yaml.safe_load(f)

work_area = card['campaign']['workArea']
if os.path.isdir(work_area):
    if raw_input("``workArea: {}`` already exists. Continue? (y/n)".format(work_area)) != "y":
        exit()
else:
    os.mkdir(work_area)

if(card['campaign']['tag_mod'] is not None) and (card['campaign']['tag_extension'] is not None):
    print("Can't specify both ``campaign: tag_mod`` and ``campaign: tag_extension``. Leave one empty.")
    exit()

with open(card['campaign']['crab_template'], 'r') as template_file:
    base_crab_config = template_file.read()

if card['campaign']['datasets'].endswith(".txt"):
    with open(card['campaign']['datasets'], 'r') as dataset_file:
        datasets = [d for d in dataset_file.read().split() if len(d) > 10 and not d.startswith("#")]
else:
    datasets = [d for d in card['campaign']['datasets'].split("\n") if len(d) > 10 and not d.startswith("#")]

if args.make:
    print("Making configs:")
    for dataset in datasets:
        print("   ==> "+dataset)
        crab_config = copy.deepcopy(base_crab_config)
        dataset_name = dataset.lstrip("/").replace("/", "_")

        tag = dataset.split("/")[2]
        if card['campaign']['tag_mod'] is not None:
            tag = re.sub(r'MiniAOD[v]?[0-9]?', card['campaign']['tag_mod'], tag) if tag.startswith('RunII') else tag + '_' + card['campaign']['tag_mod']
        elif card['campaign']['tag_extension'] is not None:
            tag = tag + '_' + card['campaign']['tag_extension']
        else:
            raise ValueError("Either ``campaign: tag_mod`` or ``campaign: tag_extension`` need to be specified")

        verbatim_lines = []
        card_info = {
            '_requestName_': dataset_name,
            '_workArea_': card['campaign']['workArea'],
            '_psetName_': card['campaign']['config'],
            '_inputDataset_': dataset,
            '_outLFNDirBase_': card['campaign']['outLFNDirBase'],
            '_storageSite_': card['campaign']['storageSite'],
            '_publication_': str(card['campaign']['publication']),
            '_splitting_': 'LumiBased' if card['campaign']['data'] else "Automatic",
            '_outputDatasetTag_': tag, 
        }
        
        if args.test:
            verbatim_lines.append("config.Data.totalUnits = 1")
            card_info['_publication_'] = 'False'

        if card['campaign']['data']:
            verbatim_lines.append("config.Data.unitsPerJob = 50")
            verbatim_lines.append("config.JobType.maxJobRuntimeMin = 2750")
        if card['campaign']['data'] and card['campaign']['lumiMask'] is not None:
            verbatim_lines.append("config.Data.lumiMask = '{}'".format(card['campaign']['lumiMask']))
        if card['campaign']['voGroup'] is not None:
            verbatim_lines.append("config.User.voGroup = '{}'".format(card['campaign']['voGroup']))
        
        for line in verbatim_lines:
            crab_config += "\n" + line
        crab_config += "\n"

        for key in card_info:
            crab_config = crab_config.replace(key, card_info[key])

        cfg_filename = os.path.join(card['campaign']['workArea'] , 'submit_{}.py'.format(dataset_name))
        with open(cfg_filename, 'w') as cfg_file:
            cfg_file.write(crab_config)


if args.submit:
    from multiprocessing import Process
    import importlib
    print("Submitting configs:")
    for dataset in datasets:
        print("   ==> "+dataset)
        dataset_name = dataset.lstrip("/").replace("/", "_")
        cfg_filename = os.path.join(card['campaign']['workArea'] , 'submit_{}.py'.format(dataset_name))
        config_base = importlib.import_module(cfg_filename)
        p = Process(target=submit, args=(config_base.config,))
        p.start()
        p.join()
