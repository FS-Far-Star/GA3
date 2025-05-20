def b1(Re):
    if Re < 10:
        return 48.0
    elif Re > 10 and Re < 100:
        return 45.1
    elif Re > 100 and Re < 1000:
        return 4.57
    elif Re > 1000 and Re < 10000:
        return 0.468
    elif Re > 10000 and Re < 100000:
        return 0.372
    
def b2(Re):
    if Re < 10:
        return -1.0
    elif Re > 10 and Re < 100:
        return -0.973
    elif Re > 100 and Re < 1000:
        return -0.476
    elif Re > 1000 and Re < 10000:
        return -0.152
    elif Re > 10000 and Re < 100000:
        return -0.123
    
def b3(Re):
    return 7.0

def b4(Re):
    return 0.5