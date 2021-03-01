import os

# -------------------- reversed version ------------------- #
# Sources : 
# https://gist.github.com/jrivero/1085501
# https://stackoverflow.com/questions/10933838/how-to-read-a-csv-file-in-reverse-order-in-python

def reversed_lines(file):
    "Generate the lines of file in reverse order."
    part = ''
    for block in reversed_blocks(file):
        for c in reversed(block):
            if c == '\n' and part:
                yield part[::-1]
                part = ''
            part += c
    if part: yield part[::-1]

def reversed_blocks(file, blocksize=4096):
    "Generate blocks of file's contents in reverse order."
    file.seek(0, os.SEEK_END)
    here = file.tell()
    while 0 < here:
        delta = min(blocksize, here)
        here -= delta
        file.seek(here, os.SEEK_SET)
        yield file.read(delta)

def split(filehandler, delimiter=',', iterations = 500, index = -1, reverse = True, 
    output_name_template='output_%s.csv', output_path='.', keep_headers=True, expected_lenght = 9, idx_to_delete = [5,6]):
    """
    Splits a CSV file into multiple pieces.
    
    A quick bastardization of the Python CSV library.
    Arguments:
        `row_limit`: The number of rows you want in each output file. 10,000 by default.
        `output_name_template`: A %s-style template for the numbered output files.
        `output_path`: Where to stick the output files.
        `keep_headers`: Whether or not to print the headers in each output file.
    Example usage:
        >> from toolbox import csv_splitter;
        >> csv_splitter.split(open('/home/ben/input.csv', 'r'));
    """
    import csv
    if(reverse):
        reader = csv.reader(reversed_lines(filehandler), delimiter=delimiter)
    else:
        reader = csv.reader(filehandler, delimiter=delimiter)
    current_piece = 1
    current_out_path = os.path.join(
         output_path,
         output_name_template  % current_piece
    )
    current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
    current_limit = iterations

    if reverse:
        if keep_headers:
            headers = next(csv.reader(filehandler, delimiter=delimiter))
        else:
            headers = next(reader)
        current_out_writer.writerow(headers)
    max_it = -1
    l = [0,0,0,0,0,0]
    for i, row in enumerate(reader):
        if(len(row)<expected_lenght):
            l[0]+=1
            continue
        elif(len(row)==expected_lenght+1):
            l[1]+=1
            row = [row[r] for r in range(len(row)) if r not in [5]]
        elif(len(row)==expected_lenght+2):
            l[2]+=1
            row = [row[r] for r in range(len(row)) if r not in [5,6]]
        elif(len(row)==expected_lenght+3):
            l[3]+=1
            row = [row[r] for r in range(len(row)) if r not in [1,2,3]]
        elif(len(row)!=expected_lenght):
            l[4]+=1
            print(row)
        else :
            l[5]+=1
        try :
            it = int(row[index]) # idx of liens should be given
        except ValueError:
            continue
        if(i==0):
            max_it = it
            print(max_it)
        if it < max_it - current_limit:
            current_piece += 1
            current_limit = iterations * current_piece
            current_out_path = os.path.join(
               output_path,
               output_name_template  % current_piece
            )
            current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
            if keep_headers:
                current_out_writer.writerow(headers)
        current_out_writer.writerow(row)
    print('Under / +1 / +2 / +3 / else / normal : {}'.format(l))