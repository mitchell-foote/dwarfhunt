p_l = 2 

p_w = 3

p_h = 4

def total_square_feet( l, w, h): 
    return 2*l*w + 2*w*h + 2*h*l + min(l, w, h)

print(total_square_feet(p_l, p_w, p_h))