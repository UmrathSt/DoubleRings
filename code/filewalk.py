import os

def fileList(source, beginswith, endswith):
    """Return a list of filenames locted in the folder SOURCE
       which begin with BEGINSWITH and end with the string
       endswith
    """
    matches = []
    for root, dirnames, filenames in os.walk(source, topdown=False):
        for filename in filenames:
            if filename.endswith((endswith)) and filename.startswith((beginswith)):
                matches.append(os.path.join(root, filename))
    return sorted(matches)

