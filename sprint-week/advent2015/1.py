floor = 0

targetString = "((())"

for char in targetString:
    if char == "(":
        floor += 1
    elif char == ")":
        floor -= 1

print(floor)

