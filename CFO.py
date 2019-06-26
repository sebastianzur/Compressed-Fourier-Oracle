from projectq import MainEngine
from projectq.ops import X, C, Measure, All, Swap
from projectq.meta import Control
import random
import sys

# Locate X query in the database
def Locate(eng, X_reg, D_reg, L_reg):
    
    print("Running Locate")
    
    # Get registers of database
    D_reg_X = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) < m)]
    D_reg_Y = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) >= m)]
    
    # Initialize auxiliary registers A,B
    A_reg = eng.allocate_qureg(m)
    B_reg = eng.allocate_qubit()
    
    All(X) | A_reg
    All(X) | B_reg
    
    for i in range(q):
        
        # Get registers of database
        D_reg_X_i = [qubit for qubit in D_reg_X if (D_reg_X.index(qubit) >= m * i \
                     and D_reg_X.index(qubit) < m * (i + 1))]
        D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                     and D_reg_Y.index(qubit) < n * (i + 1))]
        
        # Compute difference
        for j in range(m):
            C(X, 1) | (X_reg[j], A_reg[j])
            C(X, 1) | (D_reg_X_i[j], A_reg[j])
            
        # Check if D_i^Y =/= 0 and save location if difference is also 0
        All(X) | D_reg_Y_i
        C(X, n) | (D_reg_Y_i, B_reg)
        C(X, m + 1) | (B_reg + A_reg, L_reg[i])
        C(X, n) | (D_reg_Y_i, B_reg)
        All(X) | D_reg_Y_i
        
        # Uncompute A
        for j in range(m):
            C(X, 1) | (X_reg[j], A_reg[j])
            C(X, 1) | (D_reg_X_i[j], A_reg[j])
    
    # Clean up
    All(X) | A_reg
    All(X) | B_reg
    del A_reg
    del B_reg
    print("Finished Locate")        
    
    return L_reg

# Check if the qubits saved in register X are bitwise larger than those in register Y
# Output the Boolean result in register R
def Larger(eng, U_reg, V_reg, R_reg):
    
    # Initialize auxiliary register A
    t = len(U_reg)
    A_reg = eng.allocate_qureg(3 * t - 1)
    
    for i in range(t - 1):
        
        # Save if X_i > Y_i
        X | V_reg[i]
        C(X, 2) | (U_reg[i], V_reg[i], A_reg[3 * i])
        X | V_reg[i]
        
        # Save if Y_i > X_i
        X | U_reg[i]
        C(X, 2) | (U_reg[i], V_reg[i], A_reg[3 * i + 1])
        X | U_reg[i]
        
        # Save if X_i = Y_i
        X | A_reg[3 * i]
        X | A_reg[3 * i + 1]
        C(X, 2) | (A_reg[3 * i], A_reg[3 * i + 1], A_reg[3 * i + 2])
        X | A_reg[3 * i]
        X | A_reg[3 * i + 1] 
    
    # Save if X_t > Y_t    
    X | V_reg[t - 1]
    C(X, 2) | (U_reg[t - 1], V_reg[t - 1], A_reg[3 * (t - 1)])
    X | V_reg[t - 1]
    
    # Save if Y_t > X_t
    X | U_reg[t - 1]
    C(X, 2) | (U_reg[t - 1], V_reg[t - 1], A_reg[3 * t - 2])
    X | U_reg[t - 1]
    
    # Locate highest bit difference
    for i in reversed(range(t - 1)):
        C(X, 2) | (A_reg[3 * i + 2], A_reg[3 * (i + 1)], A_reg[3 * i])
        C(X, 2) | (A_reg[3 * i + 2], A_reg[3 * (i + 1) + 1], A_reg[3 * i + 1])
    
    # Save result    
    C(X, 1) | (A_reg[0], R_reg)
    
    # Uncompute everything by doing the calculation in reverse order
    for i in reversed(range(t - 1)):
        C(X, 2) | (A_reg[3 * i + 2], A_reg[3 * (i + 1)], A_reg[3 * i])
        C(X, 2) | (A_reg[3 * i + 2], A_reg[3 * (i + 1) + 1], A_reg[3 * i + 1])
        
    X | V_reg[t - 1]
    C(X, 2) | (U_reg[t - 1], V_reg[t - 1], A_reg[3 * (t - 1)])
    X | V_reg[t - 1]
    
    X | U_reg[t - 1]
    C(X, 2) | (U_reg[t - 1], V_reg[t - 1], A_reg[3 * t - 2])
    X | U_reg[t - 1]
    
    for i in reversed(range(t - 1)):
        X | V_reg[i]
        C(X, 2) | (U_reg[i], V_reg[i], A_reg[3 * i])
        X | V_reg[i]
        
        X | U_reg[i]
        C(X, 2) | (U_reg[i], V_reg[i], A_reg[3 * i + 1])
        X | U_reg[i]
        
        X | A_reg[3 * i]
        X | A_reg[3 * i + 1]
        C(X, 2) | (A_reg[3 * i], A_reg[3 * i + 1], A_reg[3 * i + 2])
        X | A_reg[3 * i]
        X | A_reg[3 * i + 1] 
    
    # Clean up
    del A_reg
    
    return R_reg

# Permute the database into the correct order, based on register A
def Permute(eng, D_reg, A_reg):
    
    print("Running Permute")
    for i in range(q):
        with Control(eng, A_reg[i]):
            for j in range(i, q - 1):
                for k in range(m + n):
                    Swap | (D_reg[j * (m + n) + k], D_reg[(j + 1) * (m + n) + k])
                    
    print("Finished Permute")
    
    return D_reg

# Inverse of above permutation
def Permute_inv(eng, D_reg, A_reg):
    
    print("Starting Permute_inv")
    for i in reversed(range(q)):
        with Control(eng, A_reg[i]):
            for j in reversed(range(i, q - 1)):
                for k in reversed(range(m + n)):
                    Swap | (D_reg[(j + 1) * (m + n) + k], D_reg[j * (m + n) + k])
    
    print("Finished Permute_inv")
                
    return D_reg

# Add query to the database
def Add(eng, X_reg, D_reg, L_reg):
    
    print("Running Add")
    
    # Initialize auxiliary register A
    A_reg = eng.allocate_qureg(q)
    
    # Get registers of database
    D_reg_X = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) < m)]
    D_reg_Y = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) >= m)]
    
    for i in range(q):
        
        # Get registers of database
        D_reg_X_i = [qubit for qubit in D_reg_X if (D_reg_X.index(qubit) >= m * i \
                     and D_reg_X.index(qubit) < m * (i + 1))]
        D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                     and D_reg_Y.index(qubit) < n * (i + 1))]
        
        # Look for where the database entries are larger than the query
        A_reg[i] = Larger(eng, D_reg_X_i, X_reg, A_reg[i])
        
        # Correct for empty entries
        All(X) | D_reg_Y_i
        C(X, n) | (D_reg_Y_i, A_reg[i])
        All(X) | D_reg_Y_i
        
        # Flip all higher bits, such that Hamming weight of A becomes 1
        for j in range(i + 1, q):
            C(X, 1) | (A_reg[i], A_reg[j])
            
    # Permute the database into the correct order
    D_reg = Permute_inv(eng, D_reg, A_reg)
    
    for i in range(q):
        with Control(eng, A_reg[i]):
            
            # Get registers of database
            D_reg_X_i = [qubit for qubit in D_reg_X if (D_reg_X.index(qubit) >= m * i \
                         and D_reg_X.index(qubit) < m * (i + 1))]
            
            # Add to the database and update register L
            for j in range(m):
                C(X, 1) | (X_reg[j], D_reg_X_i[j])
            X | L_reg[i] 
    
    # Clean Up
    All(X) | X_reg
    with Control(eng, X_reg):
        for i in range(q):
            C(X, 1) | (L_reg[i], A_reg[i])
    All(X) | X_reg
    del A_reg
    print("Finished Add")            
    
    return D_reg, L_reg

# XOR the Y query into the database into the correct location
def Update(eng, Y_reg, D_reg, L_reg):
    
    print("Running Update")
    D_reg_Y = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) >= m)]
    for i in range(q):
        
        # Get registers of database
        D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                     and D_reg_Y.index(qubit) < n * (i + 1))]
        for j in range(n):
            C(X, 2) | (Y_reg[j], L_reg[i], D_reg_Y_i[j])
    print("Finished Update")
    
    return D_reg

# Remove invalid entries from the database
def Remove(eng, X_reg, D_reg, L_reg):
    
    print("Running Remove")
    
    # Initialize auxiliary register A, B
    A_reg = eng.allocate_qureg(q)
    B_reg = eng.allocate_qubit()
    
    # Get registers of database
    D_reg_X = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) < m)]
    D_reg_Y = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) >= m)]
    
    for i in range(q):
        # Get registers of database
        D_reg_X_i = [qubit for qubit in D_reg_X if (D_reg_X.index(qubit) >= m * i \
                     and D_reg_X.index(qubit) < m * (i + 1))]
        D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                     and D_reg_Y.index(qubit) < n * (i + 1))]
        
        # Remove entry if D_i^Y = 0 and entry has been modified during the query
        All(X) | D_reg_Y_i
        for j in range(m) :
            C(X, n + 2) | (D_reg_Y_i + [L_reg[i]] + [X_reg[j]], D_reg_X_i[j])
        
        # Save that we have removed an entry
        C(X, n + 1) | (D_reg_Y_i + [L_reg[i]], B_reg)
        All(X) | D_reg_Y_i
    
    # If entry has been removed
    with Control(eng, B_reg):   
        for i in range(q):
            
            # Get registers of database
            D_reg_X_i = [qubit for qubit in D_reg_X if (D_reg_X.index(qubit) >= m * i \
                         and D_reg_X.index(qubit) < m * (i + 1))]
            D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                         and D_reg_Y.index(qubit) < n * (i + 1))]
            
            # Find locations where query is larger than the database
            A_reg[i] = Larger(eng, X_reg, D_reg_X_i, A_reg[i])
            
            # Correct for empty entries which have not been modified during the query
            X | L_reg[i]
            All(X) | D_reg_Y_i
            C(X, n + 1) | (D_reg_Y_i + [L_reg[i]], A_reg[i])
            All(X) | D_reg_Y_i
            X | L_reg[i]
            
            # Flip all lower bits, such that Hamming weight of A becomes 1 
            for j in reversed(range(i)):
                C(X, 1) | (A_reg[i], A_reg[j])
       
        # Update the location register        
        for i in range(q):
            C(X, 1) | (A_reg[i], L_reg[i])
        
        # Permutate databse into the correct order
        D_reg = Permute(eng, D_reg, A_reg)
        
        # Uncompute A by doing the A computation of the Add subroutien
        for i in reversed(range(q)):
            
            # Get registers of databas3
            D_reg_X_i = [qubit for qubit in D_reg_X if (D_reg_X.index(qubit) >= m * i \
                         and D_reg_X.index(qubit) < m * (i + 1))]
            D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                         and D_reg_Y.index(qubit) < n * (i + 1))]
            
            for j in reversed(range(i, q)):
                C(X, 1) | (A_reg[i], A_reg[j])
            
            # Correct for empty entries
            All(X) | D_reg_Y_i
            C(X, n) | (D_reg_Y_i, A_reg[i])
            All(X) | D_reg_Y_i
           
            A_reg[i] = Larger(eng, D_reg_X_i, X_reg, A_reg[i])
    
    # Reset register B
    A_reg = Locate(eng, X_reg, D_reg, A_reg)
    All(X) | A_reg
    C(X, q) | (A_reg, B_reg)
    All(X) | A_reg
    A_reg = Locate(eng, X_reg, D_reg, A_reg)
    
    # Clean up
    del A_reg
    del B_reg
    
    print("Finished Remove")
    
    return D_reg, L_reg

# Clean up the auxiliary A register
def Cleanup(eng, Y_reg, D_reg, L_reg, A_reg):
    
    print("Running Locate")
    
    # Get registers of database
    D_reg_Y = [qubit for qubit in D_reg if (D_reg.index(qubit) % (m + n) >= m)]
    
    # Initialize auxiliary register B
    B_reg = eng.allocate_qureg(n)
    All(X) | B_reg
    
    for i in range(q):
        
        # Get registers of database
        D_reg_Y_i = [qubit for qubit in D_reg_Y if (D_reg_Y.index(qubit) >= n * i \
                     and D_reg_Y.index(qubit) < n * (i + 1))]
        
        # Compute difference in Y values
        for j in range(n):
            C(X, 1) | (Y_reg[j], B_reg[j])
            C(X, 1) | (D_reg_Y_i[j], B_reg[j])
            
        # Check if D_i equals the query
        C(X, n + 1) | ([L_reg[i]] + B_reg, A_reg)
        
        # Uncompute B
        for j in range(n):
            C(X, 1) | (Y_reg[j], B_reg[j])
            C(X, 1) | (D_reg_Y_i[j], B_reg[j])
    
    # Clean up
    All(X) | B_reg
    del B_reg
    
    print("Finished Locate")        
    return A_reg

# Test a random single query to the database
def Run_test(eng):
    
    max_qubits = q * (n + m + 2) + (4 * m) + n + 1
    
    # Randomly fill register X
    X_input = ''.join(str(random.randint(0, 1)) for _ in range(m))
    X_reg = eng.allocate_qureg(m)
    for i in range(m):
        if X_input[i] == '1':
            X | X_reg[i]
        
    # Randomly fill register Y
    Y_input = ''.join('0' for _ in range(n))
    zero = Y_input
    # Zero not allowed
    while Y_input == zero:
        Y_input = ''.join(str(random.randint(0, 1)) for _ in range(n))
    Y_reg = eng.allocate_qureg(n)
    for i in range(n):
        if Y_input[i] == '1':
            X | Y_reg[i]
            
    # Randomly fill register D, apart from one empty entry
    # Generate database x entries and order them
    database_x = []
    database_y = []
    for i in range(0, q - 1):
        D_input = ''.join(str(random.randint(0, 1)) for _ in range(m))
        while D_input in database_x:
            D_input = ''.join(str(random.randint(0, 1)) for _ in range(m))
        database_x.append(D_input)
        D_input = ''.join(str(random.randint(0, 1)) for _ in range(n))
        # Zero not allowed
        while D_input == zero:
            D_input = ''.join(str(random.randint(0, 1)) for _ in range(n))
        database_y.append(D_input)
    database_x.sort()
    
    # Add one empty entry
    database_x.append(''.join('0' for _ in range(m)))
    database_y.append(''.join('0' for _ in range(n)))
        
    D_reg = eng.allocate_qureg((m + n) * q)
    
    for i in range(0, q - 1):
        for j in range(m):
            if database_x[i][j] == '1':
                X | D_reg[(m + n) * i + j]
        for j in range(n):
            if database_y[i][j] == '1':
                X | D_reg[(m + n) * i + j + m]
    
    # Auxilliary registers A, L
    A_reg = eng.allocate_qubit()
    L_reg = eng.allocate_qureg(q)
    
    # Locate X in the database
    L_reg = Locate(eng, X_reg, D_reg, L_reg)
    
    # Add X to the database if not located 
    All(X) | L_reg
    C(X, q) | (L_reg, A_reg)
    All(X) | L_reg
    with Control(eng, A_reg):
        D_reg, L_reg = Add(eng, X_reg, D_reg, L_reg)
    
    # Update the Y entry in the database
    D_reg = Update(eng, Y_reg, D_reg, L_reg)
    
    # Remove invalid entries
    D_reg, L_reg = Remove(eng, X_reg, D_reg, L_reg)
    
    # Reset auxiliary register A, L
    A_reg = Cleanup(eng, Y_reg, D_reg, L_reg, A_reg) 
    L_reg = Locate(eng, X_reg, D_reg, L_reg)
    
    # Measure the qubits
    All(Measure) | D_reg
    All(Measure) | X_reg
    All(Measure) | Y_reg
    
    # Execute Measurements
    eng.flush()
    
    # Print database before the query, this is classical!
    print("\nDatabase before the query:")
    for i in range(q):
        print("Entry " + str(i + 1) + ": {} | {}".format([int(database_x[i][j]) \
              for j in range(m)], [int(database_y[i][j]) for j in range(n)]))
    
    # Print joint state after the query
    print("Query" + ": {} | {}".format([int(X_reg[i]) for i in range(len(X_reg))], \
           [int(Y_reg[i]) for i in range(len(Y_reg))]))
    print("Database after the query:")
    for i in range(q):
        print("Entry " + str(i + 1) + ": {} | {}".format([int(D_reg[j]) for j \
              in range(len(D_reg)) if j >= (m + n) * i and j < (m + n) * i + m], \
              [int(D_reg[j]) for j in range(len(D_reg)) if j >= (m + n) * i + m \
              and j < (m + n) * (i + 1)]))
    print("Number of qubits used: " + str(max_qubits))
    
# Simulate multiple database queries   
def Run_sim(eng):
    
    # Register D
    D_reg = eng.allocate_qureg((m + n) * q)
    
    # Auxilliary registers A, L
    A_reg = eng.allocate_qubit()
    L_reg = eng.allocate_qureg(q)
    
    going = 1
    max_qubits = q * (n + m + 2) + (4 * m) + n + 1
    
    while going:
        
        # Register X
        X_input = input("Input the X register bitstring: ")
        X_reg = eng.allocate_qureg(m)
        for i in range(m):
            try:
                if X_input[i] == '1':
                    X | X_reg[i]
            except:
                sys.exit("faulty input, program terminated")
        
        # Register Y
        Y_input = input("Input the Y register bitstring: ")
        Y_reg = eng.allocate_qureg(n)
        for i in range(n):
            try:
                if Y_input[i] == '1':
                    X | Y_reg[i]
            except:
                sys.exit("faulty input, program terminated")
                
        # Locate X in the database
        L_reg = Locate(eng, X_reg, D_reg, L_reg)
        
        # Add X to the database if not located 
        All(X) | L_reg
        C(X, q) | (L_reg, A_reg)
        All(X) | L_reg
        with Control(eng, A_reg):
            D_reg, L_reg = Add(eng, X_reg, D_reg, L_reg)
            
        # Update the Y entry in the database
        D_reg = Update(eng, Y_reg, D_reg, L_reg)
        
        # Remove invalid entries
        D_reg, L_reg = Remove(eng, X_reg, D_reg, L_reg)
        
        # Reset auxiliary register A, L
        A_reg = Cleanup(eng, Y_reg, D_reg, L_reg, A_reg) 
        L_reg = Locate(eng, X_reg, D_reg, L_reg)
       
        # Reset the query
        del X_reg
        Y_reg_save = eng.allocate_qureg(n)
        max_qubits = max_qubits + n
        for i in range(n):
            Swap | (Y_reg_save[i], Y_reg[i])
        del Y_reg
        going = int(input("Press (0) to stop, or (1) to continue: "))
    
    # Measure the qubits
    All(Measure) | D_reg
    
    # Execute Measurements and print results
    eng.flush()
    print("\nDatabase after the query:")
    for i in range(q):
        print("Entry " + str(i + 1) + ": {} | {}".format([int(D_reg[j]) for j \
              in range(len(D_reg)) if j >= (m + n) * i and j < (m + n) * i + m], \
              [int(D_reg[j]) for j in range(len(D_reg)) if j >= (m + n) * i + m \
              and j < (m + n) * (i + 1)]))
    print("Number of qubits used: " + str(max_qubits))
    
if __name__ == "__main__":

    # Create the compiler
    eng = MainEngine()
    
    # Choose what instance to run
    choice = True
    while choice:
        try:
            run = int(input("For random test query type 1, for full simulation type 2\n"))
        except:
            sys.exit("faulty input, program terminated")
            
        if run == 1:
            choice = False
            # Sizes are fixed for random test query
            m = 2
            n = 1
            q = 3
            Run_test(eng)
        elif run == 2:
            choice = False
            try:
                m = int(input("X size: "))
                n = int(input("Y size: "))
                q = int(input("Database size: "))
            except:
                sys.exit("faulty input, program terminated")
            Run_sim(eng)
        else:
            print("Please type either 1 or 2\n")