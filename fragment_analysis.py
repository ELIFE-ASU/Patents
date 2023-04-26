import pickle
import igraph as ig


def main():
    ### Fragment count testing
    fragment_count = pickle.load(file=open(
        "Data/AssemblyValues/Fragments/authorFragsCount_1973-03_updated.p",
        "rb"))
    fragments = pickle.load(file=open(
        "Data/AssemblyValues/Fragments/authorFrags_1973-03_updated.p", "rb"))

    print(fragment_count)
    # print(fragments)

    print(ig.summary(fragments[50]))
    print(fragments[50].vs["color"])
    print(fragments[50].es["color"])


if __name__ == "__main__":
    main()
