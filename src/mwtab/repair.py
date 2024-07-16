# -*- coding: utf-8 -*-
"""
The repair command.
"""

import re
import traceback
import sys

from . import mwtab


{
 "Duplicate keys in SSF Additional sample data, and order matters." : [],
 "Internal representations do not match." : ["partially repairable"], # This is more generalized and some can be repaired.
 "WORKBENCH header has duplicate keys." : ["partially repairable"], # If the key also has duplicate values can be repaired. For different values just add a prefix maybe?
 "Non key value element(s) in WORKBENCH header missing." : [], # This is more of an issue of having non key value elements at all and whether that's an error or not.
 # TODO change parsing so it pairs up the filename before the datatrack_ID and give the filename the key of "DATATRACK_ID:###". 
 # In other wrods, make the key the associated datatrack_ID key value pair.
 # Need t othink about this and look at it some more. Would be much easier to just ignore this stuff. Not sure how much value it has.
 # TODO think about adding a more complicated check for data where you compare the type of the data in each column, and size, etc
 # and use that information to make corrections. This could work for AN000593 (search for this in temp.py it explains the situation).
 "Duplicate sample names." : [],
 # AN000144 has the samples in SSF in a different order than what is in the header. The duplicated sample is near the middle of the list 
 # in the header, S00011727, but the missing sample is S00011706 which should be near the beginning, so it could be that all the names 
 # should be shifted so that S00011706 is in the correct spot, or it could be that S00011706 should replace one of the S00011727's. 
 # I don't think we can say, so I don't think this can be repaired.
 
 # TODO Duplicate factor names in SSF. AN000379 
 # TODO Changed the tokenizer so factors split on ' | ' instead of '|'. Make sure that files still parse okay. This was to fix an issue in AN001296 where a factor had a pipe in it.
 }

def data_block_element_fix(start_index, end_index, line, index, lines):
    """
    """
    if start_index and not end_index and '"' in line:
        tokens = []
        line_copy = line
        while '\t' in line_copy:
            tab_index = line_copy.index('\t')
            token = line_copy[0:tab_index]
            if token.count('"') % 2 == 0:
                tokens.append(token)
                line_copy = line_copy[tab_index+1:]
            else:
                line_copy = line_copy[0:tab_index] + line_copy[tab_index+1:]
        tokens.append(line_copy)
        
        lines[index] = '\t'.join([token.strip('" ') for token in tokens])
        return True
    return False


def is_index_in_range(index, low_end, high_end):
    if not low_end or not high_end:
        return False
    if index > low_end and index < high_end:
        return True
    return False

def is_in_a_metabolite_block(index, data_low, data_high, extended_low, extended_high, metab_low, metab_high, binned_low, binned_high):
    if is_index_in_range(index, data_low, data_high) or\
       is_index_in_range(index, extended_low, extended_high) or\
       is_index_in_range(index, metab_low, metab_high) or\
       is_index_in_range(index, binned_low, binned_high):
        return True

def line_should_have_two_letter_code(line, index, data_low, data_high, extended_low, extended_high, metab_low, metab_high, binned_low, binned_high):
    return (not line.startswith("#") and \
            not line.startswith("VERSION") and \
            not line.startswith("CREATED_ON") and \
            not line.startswith("SUBJECT_SAMPLE_FACTORS") and \
            not "METABOLITE_DATA" in line and \
            not "BINNED_DATA" in line and \
            not line == "" and \
            not is_in_a_metabolite_block(index, 
                                         data_low, data_high, 
                                         extended_low, extended_high, 
                                         metab_low, metab_high,
                                         binned_low, binned_high))

def line_is_numbers(line):
    tokens = line.split('\t')[1:]
    bool_list = []
    for token in tokens:
        try:
            float(token)
            bool_list.append(True)
        except Exception:
            if token == '':
                bool_list.append(True)
            else:
                bool_list.append(False)
    return all(bool_list)


list_headers = {"Factors" : {"optional": True, 
                             "all_spellings": ["Factors", "factors"],
                             "index_after_start": 3}, 
                "Samples" : {"optional": False, 
                             "all_spellings": ["Samples", "samples"],
                             "index_after_start": 2}, 
                "metabolite_name" : {"optional": False, 
                                     "all_spellings": ["metabolite_name", "metabolite name", "Metabolite_name", "Metabolite_Name", "Metabolite Name", "Compound", "compound"],
                                     "index_after_start": 2}, 
                "Bin range(ppm)" : {"optional": True, 
                                    "all_spellings": ["Bin range(ppm)", "bin range(ppm)"],
                                    "index_after_start": 2}}

# list headers including misspellings
all_list_headers = [spelling for header_dict in list_headers.values() for spelling in header_dict['all_spellings']]


def repair(lines, verbose):
    """
    """
    # with open(path, encoding="utf-8") as f:
    #     file_text = f.read()
    # file_text = file_text.replace('\n\n', '\n')
    # lines = file_text.split('\n')
    
    
    # # Find METABOLITE_DATA start and end indexes
    # metabolite_data_start_index = [i for i, line in enumerate(lines) if "METABOLITE_DATA_START" in line]
    # if metabolite_data_start_index:
    #     metabolite_data_start_index = metabolite_data_start_index[0]
    
    # metabolite_data_end_index = [i for i, line in enumerate(lines) if "METABOLITE_DATA_END" in line]
    # if metabolite_data_end_index:
    #     metabolite_data_end_index = metabolite_data_end_index[0]
    
    # # Find samples.
    # if metabolite_data_start_index and \
    #    (lines[metabolite_data_start_index+1].startswith("Samples") or \
    #    lines[metabolite_data_start_index+1].startswith("metabolite_name") or \
    #    lines[metabolite_data_start_index+1].startswith("Metabolite Name")):
    #     samples = lines[metabolite_data_start_index+1].split('\t')[1:]
    #     while samples and not samples[-1]:
    #         samples.pop()
    
    
    # # Find METABOLITES start and end indexes
    # metabolites_start_index = [i for i, line in enumerate(lines) if "METABOLITES_START" in line]
    # if metabolites_start_index:
    #     metabolites_start_index = metabolites_start_index[0]
    
    # metabolites_end_index = [i for i, line in enumerate(lines) if "METABOLITES_END" in line]
    # if metabolites_end_index:
    #     metabolites_end_index = metabolites_end_index[0]
    
    # # Find metabolite headers.
    # if metabolites_start_index and \
    #    (lines[metabolites_start_index+1].startswith("Samples") or \
    #    lines[metabolites_start_index+1].startswith("metabolite_name")):
    #     metabolite_headers = lines[metabolites_start_index+1].split('\t')[1:]
    #     while metabolite_headers and not metabolite_headers[-1]:
    #         metabolite_headers.pop()
    
    # # Find EXTENDED_METABOLITE_DATA start and end indexes
    # extended_metabolite_data_start_index = [i for i, line in enumerate(lines) if "EXTENDED" in line and "_METABOLITE_DATA_START" in line]
    # if extended_metabolite_data_start_index:
    #     extended_metabolite_data_start_index = extended_metabolite_data_start_index[0]
    
    # extended_metabolite_data_end_index = [i for i, line in enumerate(lines) if "EXTENDED" in line and "_METABOLITE_DATA_END" in line]
    # if extended_metabolite_data_end_index:
    #     extended_metabolite_data_end_index = extended_metabolite_data_end_index[0]
    
    # # Find extended metabolite headers.
    # if extended_metabolite_data_start_index and \
    #    (lines[extended_metabolite_data_start_index+1].startswith("Samples") or \
    #    lines[extended_metabolite_data_start_index+1].startswith("metabolite_name")):
    #     extended_metabolite_headers = lines[extended_metabolite_data_start_index+1].split('\t')[1:]
    #     while extended_metabolite_headers and not extended_metabolite_headers[-1]:
    #         extended_metabolite_headers.pop()
    
    # # Find BINNED_DATA start and end indexes
    # binned_data_start_index = [i for i, line in enumerate(lines) if "BINNED_DATA_START" in line]
    # if binned_data_start_index:
    #     binned_data_start_index = binned_data_start_index[0]
    
    # binned_data_end_index = [i for i, line in enumerate(lines) if "BINNED_DATA_END" in line]
    # if binned_data_end_index:
    #     binned_data_end_index = binned_data_end_index[0]
    
    # # Find binned headers.
    # if binned_data_start_index and \
    #    lines[binned_data_start_index+1].startswith("Bin range(ppm)"):
    #     binned_headers = lines[binned_data_start_index+1].split('\t')[1:]
    #     while binned_headers and not binned_headers[-1]:
    #         binned_headers.pop()
    
    
    
    
    
    indexes_to_delete = []
    indexes_to_skip = []
    lines_to_add = []
    current_section = None
    current_prefix = ""
    current_section_index = None
    metabolite_data_start_index = None
    metabolites_start_index = None
    extended_metabolite_data_start_index = None
    binned_data_start_index = None
    metabolite_data_end_index = None
    metabolites_end_index = None
    extended_metabolite_data_end_index = None
    binned_data_end_index = None
    ssf_samples = 0
    for i, line in enumerate(lines):
        if i in indexes_to_skip:
            continue
        
        # Upkeep to keep track of where we are in the file.
        if line.startswith('#'):
            current_section = line[1:].strip()
            current_prefix = mwtab.MWTabFile.prefixes.get(current_section, "")
            current_section_index = i
                
        elif "EXTENDED" not in line and "METABOLITE_DATA_START" in line:
            metabolite_data_start_index = i
        elif "METABOLITES_START" in line:
            metabolites_start_index = i
        elif "EXTENDED" in line and "_METABOLITE_DATA_START" in line:
            extended_metabolite_data_start_index = i
        elif "BINNED_DATA_START" in line:
            binned_data_start_index = i
        
        elif "EXTENDED" not in line and "METABOLITE_DATA_END" in line:
            metabolite_data_end_index = i
        elif "METABOLITES_END" in line:
            metabolites_end_index = i
        elif "EXTENDED" in line and "_METABOLITE_DATA_END" in line:
            extended_metabolite_data_end_index = i
        elif "BINNED_DATA_END" in line:
            binned_data_end_index = i
        
        if line.startswith('SUBJECT_SAMPLE_FACTORS'):
            ssf_samples += 1
        
        
        # Remove tabs from the workbench header.
        if line.startswith("#METABOLOMICS WORKBENCH") and '\t' in line:
            lines[i] = re.sub(r'\t+', ' ', line)
            continue
        
        # Fix it when the header accidently is on 2 different lines.
        if i > 0 and lines[i-1].startswith("#METABOLOMICS WORKBENCH") and '\t' not in line:
                lines[i-1] = lines[i-1] + ' ' + line.lstrip()
                indexes_to_delete.append(i)
                continue
        
        # If VERSION doesn't have a value then give it one.
        if match_re := re.match(r"VERSION(.*)", line):
            version_value = match_re.group(1).strip()
            if version_value == '':
                lines[i] = "VERSION             \t1"
                continue
            else:
                lines[i] = "VERSION             \t" + version_value
                continue
            # if '\t' not in match_re.group(1):
            #     lines[i] = "VERSION             \t" + match_re.group(1).lstrip()
            #     continue
        
        # Fix tabs in the sub_section names.
        if (match_re := re.match(r'(\s*[A-Z_]+\s*[A-Z_]*:?\s*[A-Z_]+)( +\t.*)', line)) and '\t' in match_re.group(1):
            lines[i] = "".join(match_re.group(1).split()) + match_re.group(2)
            continue
        
        # Fix tabs in #SECTIONS.
        if line.startswith('#') and '\t' in line:
            lines[i] = "".join(line.split())
            continue
        
        # Fix names for list header lines.
        if metabolite_data_start_index:
            if i - metabolite_data_start_index == 1:
                # This can't be done for samples because they have to match up to the samples in SSF.
                # samples = line.split('\t')[1:]
                # while samples and not samples[-1]:
                #     samples.pop()
                
                # header_count = {}
                # for j, header in enumerate(samples):
                #     if header in header_count:
                #         samples[j] = f"{header}_{header_count[header]}"
                #         header_count[header] += 1
                #     else:
                #         header_count[header] = 1
                
                # lines[i] = '\t'.join(["Samples"] + samples)
                # continue
                
                # num_of_metabolite_data_headers = len(line.split('\t'))
                if not line.startswith("Samples"):
                    lines[i] = re.sub(r"[^\t]*\t", "Samples\t", line, count=1)
                    continue
                continue
            
            elif i - metabolite_data_start_index == 2:
                # This can't be done for factors because they have to match up to the samples in SSF.
                # factors = line.split('\t')[1:]
                # while factors and not factors[-1]:
                #     factors.pop()
                
                # header_count = {}
                # for j, header in enumerate(factors):
                #     if header in header_count:
                #         factors[j] = f"{header}_{header_count[header]}"
                #         header_count[header] += 1
                #     else:
                #         header_count[header] = 1
                
                # lines[i] = '\t'.join(["Factors"] + factors)
                # continue
                
                if not line.startswith("Factors") and not line_is_numbers(line):
                    lines[i] = re.sub(r"[^\t]*\t", "Factors\t", line, count=1)
                    continue
                continue
        
        if metabolites_start_index:
            if i - metabolites_start_index == 1:
                metabolite_headers = line.split('\t')[1:]
                # num_of_metabolites_headers = len(metabolite_headers) + 1
                while metabolite_headers and not metabolite_headers[-1]:
                    metabolite_headers.pop()
                
                if "Metabolite" in metabolite_headers:
                    metabolite_headers = [header if header != "Metabolite" else "Metabolite Name" for header in metabolite_headers]
                
                header_count = {}
                for j, header in enumerate(metabolite_headers):
                    if header in header_count:
                        metabolite_headers[j] = f"{header}_{header_count[header]}"
                        header_count[header] += 1
                    else:
                        header_count[header] = 1
                
                lines[i] = '\t'.join(["metabolite_name"] + metabolite_headers)
                continue
                
                # if not line.startswith("metabolite_name"):
                #     lines[i] = re.sub(r"[^\t]*\t", "metabolite_name\t", line, count=1)
                #     continue
        
        if extended_metabolite_data_start_index:
            if i - extended_metabolite_data_start_index == 1:
                extended_headers = line.split('\t')[1:]
                # num_of_extended_headers = len(extended_headers) + 1
                while extended_headers and not extended_headers[-1]:
                    extended_headers.pop()
                
                if "Metabolite" in extended_headers:
                    extended_headers = [header if header != "Metabolite" else "Metabolite Name" for header in extended_headers]
                
                header_count = {}
                for j, header in enumerate(extended_headers):
                    if header in header_count:
                        extended_headers[j] = f"{header}_{header_count[header]}"
                        header_count[header] += 1
                    else:
                        header_count[header] = 1
                
                lines[i] = '\t'.join(["metabolite_name"] + extended_headers)
                continue
                
                # if not line.startswith("metabolite_name"):
                #     lines[i] = re.sub(r"[^\t]*\t", "metabolite_name\t", line, count=1)
                #     continue
        
        if binned_data_start_index:
            if i - binned_data_start_index == 1:
                # This can't be done for binned_headers because they have to match up to the samples in SSF.
                # binned_headers = line.split('\t')[1:]
                # while binned_headers and not binned_headers[-1]:
                #     binned_headers.pop()
                
                # header_count = {}
                # for j, header in enumerate(binned_headers):
                #     if header in header_count:
                #         binned_headers[j] = f"{header}_{header_count[header]}"
                #         header_count[header] += 1
                #     else:
                #         header_count[header] = 1
                
                # lines[i] = '\t'.join(["Bin range(ppm)"] + binned_headers)
                # continue
                
                # num_of_binned_headers = len(line.split('\t'))
                if not line.startswith("Bin range(ppm)"):
                    lines[i] = re.sub(r"[^\t]*\t", "Bin range(ppm)\t", line, count=1)
                    continue
                continue
                
        # Find header lines in the wrong place and delete them.
        if any(line.startswith(header) for header in all_list_headers):
            # If it's the wrong index it should have taken one of the if's above this.
            # if i - metabolite_data_start_index == 1 or\
            #    i - metabolite_data_start_index == 2 or\
            #    i - metabolites_start_index == 1 or\
            #    i - extended_metabolite_data_start_index == 1 or\
            #    i - binned_data_start_index == 1:
            #     continue
            # We assume that if this line is not where we expect it to be then it shouldn't be there.
            indexes_to_delete.append(i)
        
        # Fix the issue where tabs inside of double quotes messes things up in data blocks.
        if data_block_element_fix(metabolite_data_start_index, metabolite_data_end_index, line, i, lines):
            continue
        if data_block_element_fix(metabolites_start_index, metabolites_end_index, line, i, lines):
            continue
        if data_block_element_fix(extended_metabolite_data_start_index, extended_metabolite_data_end_index, line, i, lines):
            continue
        if data_block_element_fix(binned_data_start_index, binned_data_end_index, line, i, lines):
            continue
        
        # Fix bad prefix.
        if not line.startswith(current_prefix) and not line.startswith('#'):
            # The line could have a bad prefix or no prefix.
            # Assume that if a ':' is found close enough to the start of the line it is a bad prefix.
            # "_1 CH: ..." has been seen before, so using index 6 as the "close enough".
            if ':' in line and line.index(':') < 6:
                lines[i] = re.sub(r"[^:]*:", current_prefix, line, count=1)
            # Assume if there is no tab that the line should be added to the one above it.
            elif '\t' not in line:
                lines[i-1] = lines[i-1] + line.lstrip()
                indexes_to_delete.append(i)
            else:
                lines[i] = current_prefix + line
            continue
        
        # Remove #FACTORS lines.
        if line.strip() == "#FACTORS":
            indexes_to_delete.append(i)
            continue
        
        # Remove extra #END lines.
        if line.strip() == "#END" and any(["#END" in prev_line for prev_line in lines[i-5:i]]):
            indexes_to_delete.append(i)
            continue
        
        # Fix various problems with RESULTS_FILE lines.
        if "RESULTS_FILE" in line:
            if line.startswith('#'):
                line = line[1:]
            
            if '\t\t' in line:
                line = line.replace('\t\t', '\t')
            
            if "NMR_RESULTS_FILE" in line:
                if "NM:" not in line:
                    line = "NM:" + line
                if current_section != 'NM':
                    lines_to_add.append([line, 'NM'])
                    indexes_to_delete.append(i)
            
            if "MS_RESULTS_FILE" in line:
                if "MS:" not in line:
                    line = "MS:" + line
                if current_section != 'MS':
                    lines_to_add.append([line, 'MS'])
                    indexes_to_delete.append(i)
            
            lines[i] = line
            continue
        
        # Remove duplicate sub sections and other stuff.
        if match_re := re.match(r'((\w\w):([A-Z_]+))(.*)', line):
            # Look for duplicate sub sections.
            sub_section_indexes = [index + current_section_index for index, line in enumerate(lines[current_section_index:]) if line.startswith(match_re.group(1))]
            duplicate_line_indexes = [index + current_section_index for index, line2 in enumerate(lines[current_section_index:]) if line2 == line]
            if len(duplicate_line_indexes) > 1 and len(sub_section_indexes) != len(duplicate_line_indexes):
                indexes_to_delete.extend(duplicate_line_indexes[1:])
                indexes_to_skip.extend(duplicate_line_indexes)
                continue
            
            # Add tab in if not there. 
            if '\t' not in match_re.group(4):
                # Adding the correct spacing isn't necessary because 
                # using the MWTabFile class at the end should fix that, but I am doing it anyway.
                lines[i] = match_re.group(1) + (30 - len(match_re.group(1))) * ' ' + '\t' + match_re.group(4).lstrip()
                continue
        
        # Remove out of place sub section.
        if line.startswith("PR:INSTITUTE") and "ST:INSTITUTE" in line:
            lines[i] = line[0:line.index("ST:INSTITUTE")]
            continue
    
    
    for i, index in enumerate(sorted(indexes_to_delete)):
        del lines[index - i]
    
    for line, section_name in lines_to_add:
        try:
            section_index = [index for index, line2 in enumerate(lines) if line2 == "#" + section_name][0]
        except Exception:
            continue
        lines.insert(section_index + 1, line)


    # Try to use the MWTabFile class to get an internal representation and save it out because this will standardize spacing issues.
    # Also fixes broken subsections.
    # If that can't be done then just save out the lines as is.
    file_str = '\n'.join(lines)
    
    tabfile = mwtab.MWTabFile("", compatability_mode=True)
    try:
        tabfile._build_mwtabfile(file_str)
        return tabfile.writestr('mwtab')
    except Exception as e:
        if verbose:
            traceback.print_exception(e, file=sys.stdout)
            print()
        return file_str


