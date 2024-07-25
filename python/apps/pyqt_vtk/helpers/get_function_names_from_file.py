# -*- coding: utf-8 -*-



def get_function_names_from_file(filename:str) -> list:
    import ast
    with open(filename, "r") as file:
        root = ast.parse(file.read(), filename)
    
    function_names = [node.name for node in ast.walk(root) if isinstance(node, ast.FunctionDef)]
    return function_names


if __name__ == "__main__":
    L = get_function_names_from_file('C:/Users/andre/Desktop/tiger/geometries/planar.py')
    print(L)
    