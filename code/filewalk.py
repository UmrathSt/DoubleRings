import os

def fileList(source, beginswith):
    """Return a list of filenames locted in the folder SOURCE
       which begin with BEGINSWITH and end with .txt
    """
    matches = []
    for root, dirnames, filenames in os.walk(source):
        for filename in filenames:
            if filename.endswith(('.txt')) and filename.startswith((beginswith)):
                matches.append(os.path.join(root, filename))
    return sorted(matches)

