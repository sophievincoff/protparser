def print_attributes_and_methods(obj):
    """
    Prints all the attributes and methods of this class, to help users know what's available to them. 
    """
    attributes_and_methods = dir(obj)

    # Separate attributes and methods
    attributes = [item for item in attributes_and_methods if not callable(getattr(obj, item)) and
                    (not item.startswith('__') and not item.startswith('_'))]
    methods = [item for item in attributes_and_methods if callable(getattr(obj, item)) and
                (not item.startswith('__') and not item.startswith('_'))]

    # Print attributes and methods
    print("Attributes:")
    for attribute in attributes:
        print(f"\t{attribute}")

    print("\nMethods:")
    for method in methods:
        print(f"\t{method}")