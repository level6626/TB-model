from TB.model import TB


if __name__ == "__main__":
    for i1 in range(5):
        for i2 in range(1, 5):
            for i3 in range(3):
                for i4 in range(5):
                    for i5 in range(5):
                        A = TB(i_alpha=i1, i_Tr=i2, i_Tm=i3, i_Ta=i4, i_Mr=i5)
                        A.run_model()
