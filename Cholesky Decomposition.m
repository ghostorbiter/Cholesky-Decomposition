%%...Sample for checking...%%
%%..|

    %%...My A matrix is taken from Excercise set D (Excercise D1)...\
    %%...I decided to take 4x4 tridiagonal matrix A...\
    %%...A[1,-1,0,0; -1,2,-1,0; 0,-1,2,-1; 0,0,-1,2]...\
    %%...In 3xn matrix it looks as follows...%%
    
    A = [-1,-1,-1,0; 1,2,2,2; 0,-1,-1,-1];
    
    %%...Example B is randomly taken...%%
    
    B = [21,69,34,22; 13,70,45,88];
    
    %%...Usage of function that finds Cholesky decomposition...%%
    
    [LN, LT] = Cholesky(A, length(A));
    
    %%...Checking if Cholesky decomposition indeed works...\
    %%...Matrix L (in this case LN) should equal to...\
    %%...L = [1,1,1,1; 0,-1,-1,-1]...\
    %%...And L' (in this case LT) should equal to...\
    %%...L' = [-1,-1,-1,0; 1,1,1,1]...%%
    
    disp('Simple case #0');
    disp('Cholesky decomposition is:');
    LN, LT
    
    %%...Usage of Solve function that solves XA = B...%%
    %%...Checking if Solve function indeed works...\
    %%...The Result should be 2xn matrix X (in this case SOLVED)...\
    %%...The Result should equal to...\
    %%...X = [381,360,270,146; 440,427,344,216]...%%
    
    
    disp('Solution is:');
    SOLVED = Solve(LN, LT, B, size(B,1), size(B,2))
    
    
    %%...Checking result with build in functions...%%
    
    ACheck = [1,-1,0,0; -1,2,-1,0; 0,-1,2,-1; 0,0,-1,2];
    BCheck = [21,69,34,22; 13,70,45,88];
    
    disp('Cholesky decomposition with build-in functions is:');
    chol(ACheck)
    
    disp('Solution with build-in functions is:');
    SOLVEDCheck = BCheck/ACheck
    
%%..|

%%...More complex checks...%%
%%..|
    
    disp('Other cases:');

    %%...Declaring variables that can be changed to modify the values...%%
    %%...Length is the biggest dimension of triangular matrix...\
    %%...DimensionB is the number of rows of matrix B...\
    %%...NumOfIterations is number of test cases...%%
    %%...Range is the range of entries of randomly created matrices...%%

    Range = 100;
    Length = 10;
    DimensionB = 2;
    NumOfIterations = 5;

    %%...Initializing matrices M and B of length Length with zeroes...%%
    
    M = zeros(3,Length);
    B = zeros(DimensionB,Length);
    
    %%...Setting loop counter to 1...%%
    k = 1;

    %%...Iterating and creating test cases...%%
    
    while k <= NumOfIterations
        
        %%...Modifying matrix M by inputing random integer values...%%

        for j=1:Length
            if j == Length
                M(2,j) = randi([1,Range],1,1);
                M(1,j) = randi([-Range,Range],1,1);
            else
                M(2,j) = randi([1,Range],1,1);
                M(1,j) = randi([-Range,Range],1,1);
                M(3,j+1) = M(1,j);
            end
        end
        
        %%...Modifying last entry of upper diagonal and...\
        %%...First entry of lower diagonal to 0...%%
        
        M(1,Length) = 0;
        M(3,1) = 0;
        
        %%...Modifying matrix B by inputing random integer values...%%
        
        for i=1:DimensionB
            for j=1:Length
                B(i,j) = randi([-Range,Range],1,1);
            end
        end
        
        %%...Initializing matrix MamaM with main diagonal of M...%%
        
        MamaM = diag(M(2,:));
        
        %%...Modifying matrix MamaM to input upper and lower diagonals...\
        %%...From matrix M...%%
        
        for j=1:Length
            if j == 1
                MamaM(j,j+1) = M(1,j);
            elseif j == Length
                MamaM(j,j-1) = M(3,j);
            else
                MamaM(j,j-1) = M(3,j);
                MamaM(j,j+1) = M(1,j);
            end
        end
        
        %%...Since the randomly created matrices might not be positive...\
        %%...I decided to check if they are indeed positive definite...\
        %%...Using build-in function chol for Cholesky decomposition...\
        %%...Which return an error if the provided matrix is not...\
        %%...Positive definite. In such cases I decrement the loop...\
        %%...Counter by 1 to cover the exact number of cases...%%
        
        try 
            chol(MamaM);
        catch ME
            continue;
        end
        
        %%...Outputing the final result of all 3 matrices...%%
        
        fprintf('Case # %d\n',k);
        disp('Matrices:');
        M, B, MamaM
        
        %%...Finaly solving and outputing the result...%%
        
        [LN, LT] = Cholesky(M, length(M));
        SOLVED = Solve(LN, LT, B, size(B,1), size(B,2));
        SOLVEDCheck = B/MamaM;
        
        disp('The results are:');
        SOLVED, SOLVEDCheck
        
        %%...Incrementing loop counter...%%
        k = k+1;
        
    end
    
%%..|

%%..|

    %%...Cholesky function to find lower and upper matrices of Cholesky...\
    %%...Decomposition and it takes O(n) time...%%
    %%...Function takes 2 parameters A and n and returns tuple...%%
    %%...A is a matrix for which we want to find Cholesky decomposition..%%
    %%...n is the maximum dimension of A...%%
    %%...The returning tuple is the lower and upper triangular matrices...\
    %%...That represent Cholesky decomposition...\
    %%...They are represented as 3xn matrices...%%
    
    function [MN1, MT1]=Cholesky(A, n)
    
        %%...Here I initialize 2 matrices with zeros...\
        %%...Further I will fill them with proper values...%%
        
        MN = zeros(2,n);
        MT = zeros(2,n);

        %%...Looping through the maximum dimension of A...%%
        
        for k=1:n
            
            %%...In case loop counter is equal to 1...\
            %%...I take square roots of given matrix at position 1,1...\
            %%...Which in my 3xn case is 2,1...%%
            %%...In case loop counter is NOT equal to 1...\
            %%...I take square roots of given matrix at position k,k...\
            %%...Which in my case is 2,k with substraction of squared...\
            %%...Entry in my matrix MT at position 1, k-1 ...\
            %%...Which in Cholesky decomposition matrix L^T will be...\
            %%...Positioned at k-1, k...%%
            
            %%...It is enough to take only one entry at position k-1,k...\
            %%...Because all the entries except main 3 diagonals are 0...\
            %%...So taking them in to the formula would not result...\
            %%...In anything...%%
            
            if k == 1
                MT(2,k) = sqrt(A(2,k));
            else
                MT(2,k) = sqrt(A(2,k) - MT(1, k-1).^2);
            end
            
            %%...At this part I calculate entries of my matrix at...\
            %%...Position 1,k, which in original Cholesky decomposition...\
            %%...Matrix L^T would be at position k-1,k...\
            %%...Again, taking only 1 entry from A and MT is enough...\
            %%...Because all other entries are 0...%%
            
            MT(1,k) = A(1,k)/MT(2,k);
        end
            
        %%...Looping through maximum dimension of A...%%
        %%...Here I find MN matrix which is L matrix of Cholesky...\
        %%...Decomposition...%%
        %%...It can be easily done by inverting the MT matrix that I...\
        %%...Have found earlier, but since we deal with 3xn matrix...\
        %%...We just have to loop and assign backwards the entries of...\
        %%...MT matrix -> 
            %%...For example entry at position 2,j of MT matrix...\
            %%...Which represents main diagonal should be at position...\
            %%...1,j of MN matrix since transposing main diagonal gives...\
            %%...Same diagonal...%%
            %%...However, entries at positon 1,j of MT matrix...\
            %%...Which represents diagonal above main should be at...\
            %%...Position 2,j+1 of MN matrix since transposing...\
            %%...Diagonals except main diagonal is same as putting...\
            %%...MT's upper diagonal to MN's lower diagonal...%%
            
        for j=1:n
            
            if j == n
                MN(1,j) = MT(2,j);
            else
                MN(1,j) = MT(2,j);
                MN(2,j+1) = MT(1,j);
            end 
            
        end
        
        %%...Assigning output values to the one I find...%%
        
        MN1 = MN;
        MT1 = MT;
        
    end
   
%%..|

%%|_______________________________________|%%
%%|..........Solution for XA = B..........|%%
%%|.......................................|%%
%%|............Works for AX = B...........|%%
%%|...but need to change for loop a bit...|%%
%%|_______________________________________|%%


%%..|
    
    %%...Solve function to find the XA = B solution...\
    %%...Which takes O(n*m) time, where n and m are dimensions of B...%%
    %%...Function takes L and L^T matrices as MN and MR respectevly...\
    %%...Function also takes matrix B and its dimensions...%%
    %%...Function outputs a solution matrix, in our case matrix X...%%

    function array = Solve(MN, MT, B, m, n)
    
        %%...Here I initialize the matrix y with zeros...\
        %%...We need it to find XLL^T = B, where we replace XL^T = y...\
        %%...Now we need to solve yL = B...%%
        
        y = zeros(n,m);

        %%...Looping through dimensions of B...%%
        
        for j=1:m
            for i=1:n
                
                %%...In case loop counter is equal to 1...\
                %%...I divide the entry of B at position j,i by...\
                %%...The entry of L matrix at position i,i (in our case...\
                %%...We take MN matrix at position 1,i)...%%
                
                %%...The math in this part is easy, because all we have...\
                %%...To do is assign the appropriate values...%%
                %%...Since we have lower triangular tridiagonal matrix...\
                %%...We need to subtract previous main diagonal entry...\
                %%...From B with multiplying it by value at position...\
                %%...j,i-1 which is 2,i in our case and divide it by...\
                %%...Value at position i,i which is 1,i in our case...%%
                
                if i==1
                    y(i,j) = B(j,i)/MN(1,i);
                else
                    y(i,j) = (B(j,i)-MN(2,i)*y(i-1,j))/MN(1,i);
                end
                
            end
        end

        %%...Here I initialize the matrix X with zeros...\
        %%...We need it to find XL^T = y...%%
        
        x = zeros(m,n);
        
        %%...Looping through dimensions of B, but from bottom to top...%%
        for j=1:m
            for i=n:-1:1
                
                %%...In case loop counter is equal to n...\
                %%...I divide the entry of y at position j,i by...\
                %%...The entry of L^T matrix at position i,i (in our...\
                %%...Case we take MT matrix at position 2,i)...%%
                
                %%...The math in this part is easy as well, because...\
                %%...We need to do the same as in previous nested loop...\
                %%...Except we do it backwards...%%
                %%...Since we have upper triangular tridiagonal matrix...\
                %%...We need to substract next entry of main diagonal...\
                %%...From B with multiplying it by value at position...\
                %%...j,i+1 which is 1,i in our case and divide it by...\
                %%...Value at position i,i which is 2,i in our case...%%
                
                if i==n
                    x(j,i) = y(i,j)/MT(2,i);
                else
                    x(j,i) = (y(i,j)-MT(1,i)*x(j,i+1))/MT(2,i);
                end
                
            end
        end

        %%...Here I assign the output with x that we have found...%%
        
        array = x;
    end
    
%%..|
