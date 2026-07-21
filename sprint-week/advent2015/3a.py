
def santa_map(path):
    x = 0
    y = 0
    arivedMap = set({(0,0)})
    for char in path:
        if char == '>':
            x += 1
        elif char == '^':
            y += 1
        elif char == '<':
            x -= 1
        else:
            y -= 1
        arivedMap.add((x,y))
    return arivedMap.__len__()

print(santa_map('>'))
print(santa_map("^>v<"))
print(santa_map('^v^v^v^v^v'))