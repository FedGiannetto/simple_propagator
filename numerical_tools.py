import math

## Just a function to make things easier

def norm(vec):
    if len(vec) == 2:
        mag = math.sqrt(vec[0]**2 + vec[1]**2)
    elif len(vec) == 3:
        mag = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

    return mag