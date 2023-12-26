import requests
import sys
"""Downloads the raw unnormalized counts for a given study id
    Raises:
        Exception: if a connection issue occurs
    Returns:
        DataFrame: counts matrix
"""
def downloadCounts(study_id, fname):
        # the baseline query for remote url requests
        data_query = 'https://genelab-data.ndc.nasa.gov/genelab/data/glds/files/{}'.format(study_id)
        resp = requests.get(data_query)
        # error handling for genelab api request
        if resp.status_code != 200:
            raise Exception('GET /genelab/data {}'.format(resp.status_code))
        # filter the study files for counts
        study_response = resp.json()
        filt_file = filter(lambda x: 'Unnormalized' in x['file_name'], study_response['studies']['OSD-' + study_id]['study_files'])
        counts_remote = list(filt_file)[0]
        counts_fname = counts_remote['file_name']
        # download the counts file from genelab
        response = requests.get("https://genelab-data.ndc.nasa.gov" + counts_remote['remote_url'], stream=True)
        total = response.headers.get('content-length') # for progress bar
        with open(fname, 'wb') as file:
            total = int(total)
            downloaded = 0
            for data in response.iter_content(chunk_size=max(int(total/1000), 1024*1024)):
                downloaded += len(data)
                file.write(data)
                done = int(50*downloaded/total)
                sys.stdout.write('\r[{}{}]'.format('█' * done, '.' * (50-done)))
                sys.stdout.flush()
        return study_response
    
def downloadCountsMetadata(study_id, fname):
        # the baseline query for remote url requests
        data_query = 'https://genelab-data.ndc.nasa.gov/genelab/data/glds/files/{}'.format(study_id)
        resp = requests.get(data_query)
        # error handling for genelab api request
        if resp.status_code != 200:
            raise Exception('GET /genelab/data {}'.format(resp.status_code))
        # filter the study files for counts
        study_response = resp.json()
        filt_file = filter(lambda x: 'Metdata' in x['category'], study_response['studies']['OSD-' + study_id]['study_files'])
        counts_remote = list(filt_file)[0]
        counts_fname = counts_remote['file_name']
        # download the counts file from genelab
        response = requests.get("https://genelab-data.ndc.nasa.gov" + counts_remote['remote_url'], stream=True)
        total = response.headers.get('content-length') # for progress bar
        with open(fname, 'wb') as file:
            total = int(total)
            downloaded = 0
            for data in response.iter_content(chunk_size=max(int(total/1000), 1024*1024)):
                downloaded += len(data)
                file.write(data)
                done = int(50*downloaded/total)
                sys.stdout.write('\r[{}{}]'.format('█' * done, '.' * (50-done)))
                sys.stdout.flush()
        return study_response
    
    
