using System;
using System.Collections.Generic;
using System.Text;

// Constructor
// Copy
// RowValues
// RowVector
// ColumnValues
// ColumnVector
// SwapRows
// SwapColumns
// AddRow
// RemoveRow
// AddColumn
// RemoveColumn
// Convert1DArrTo2D
// AddArrays
// MultiplyArray
// DivideArray
// Display
// Identity (I)
// Transpose (T)
// +
// -
// DotProduct
// *
// ^ : need to support fractional powers
// ++
// --
// ReducedEchelonForm (REF)
// RowReducedEchelonForm (RREF)
// PLUDecomposition (PLU)
// IN PROGRESS:
// Inverse

namespace Matrix
{
    class Matrix
    {
        private int rows;
        private int columns;
        private double[,] values;

        public Matrix(double[,] matrixValues)
        {
            rows = matrixValues.GetLength(0);
            columns = matrixValues.GetLength(1);
            values = matrixValues;
        }

        public Matrix(int row, int col)
        {
            rows = row;
            columns = col;
            double[,] value = new double[rows, columns];
            values = value;
        }

        public Matrix()
        {
            rows = 0;
            columns = 0;
        }

        public static Matrix Copy (Matrix matrix)
        {
            Matrix copy = new Matrix(matrix.rows, matrix.columns);
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    copy.values[rowIndex, colIndex] = matrix.values[rowIndex, colIndex];
                }
            }
            return copy;
        }

        public double[] RowValues (uint index)
        {
            if (index >= rows)
            {
                Console.WriteLine("ERROR: RVal Index out of bounds...");
                return null;
            }
            double[] rowVector = new double[columns];
            for (uint colIndex = 0; colIndex < columns; colIndex++)
            {
                rowVector[colIndex] = values[index, colIndex];
            }
            return rowVector;
        }

        public Matrix RowVector(uint index)
        {
            if (index >= rows)
            {
                Console.WriteLine("ERROR: RVec Index out of bounds...");
                return null;
            }

            double[,] rowVector = new double[1, columns];
            for (uint colIndex = 0; colIndex < columns; colIndex++)
            {
                rowVector[0,colIndex] = values[index, colIndex];
            }
            Matrix vector = new Matrix(rowVector);
            return vector;
        }

        public double[] ColumnValues (uint index)
        {
            if (index >= columns)
            {
                Console.WriteLine("ERROR: CVal Index out of bounds...");
                return null;
            }
            double[] colVector = new double[rows];
            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                colVector[rowIndex] = values[rowIndex, index];
            }
            return colVector;
        }

        public Matrix ColumnVector(uint index)
        {
            if (index >= columns)
            {
                Console.WriteLine("ERROR: CVec Index out of bounds...");
                return null;
            }

            double[,] colVector = new double[rows, 1];
            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                colVector[rowIndex,0] = values[rowIndex,index];
            }
            Matrix vector = new Matrix(colVector);
            return vector;
        }

        public void SwapRows(uint row1Index, uint row2Index)
        {
            if (row1Index >= rows || row2Index >= rows)
            {
                Console.WriteLine("ERROR: SR1 Index out of bounds...");
                return;
            }

            double temp = 0;
            for (uint colIndex = 0; colIndex < columns; colIndex++)
            {
                temp = values[row1Index, colIndex];
                values[row1Index, colIndex] = values[row2Index, colIndex];
                values[row2Index, colIndex] = temp;
            }
        }

        public void SwapRows(uint rowIndex, double[] newRow)
        {
            if (rowIndex >= rows)
            {
                Console.WriteLine("ERROR: SR2 Index out of bounds...");
                return;
            }

            int vectorLength = newRow.Length;
            if (columns == 0)
            {
                columns = vectorLength;
            }

            if (vectorLength > columns)
            {
                Console.WriteLine("Warning: Input row has too many elements; truncating vector...");
            }
            else if (vectorLength < columns)
            {
                Console.WriteLine("Warning: Input row has too few elements; missing elements will be substituted with 0...");
                double[] subRow = new double[columns];
                for (uint colIndex = 0; colIndex < columns; colIndex++)
                {
                    if (colIndex < vectorLength)
                    {
                        subRow[colIndex] = newRow[colIndex];
                    }
                    else
                    {
                        subRow[colIndex] = 0;
                    }
                }
                newRow = subRow;
            }

            for (uint colIndex = 0; colIndex < columns; colIndex++)
            {
                values[rowIndex, colIndex] = newRow[colIndex];
            }
        }

        public void SwapColumns(uint col1Index, uint col2Index)
        {
            if (col1Index >= columns || col2Index >= columns)
            {
                Console.WriteLine("ERROR: SC Index out of bounds...");
                return;
            }

            double temp = 0;
            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                temp = values[rowIndex, col1Index];
                values[rowIndex, col1Index] = values[rowIndex, col2Index];
                values[rowIndex, col2Index] = temp;
            }
        }

        public void SwapColumns(uint colIndex, double[] newColumn)
        {
            if (colIndex >= columns)
            {
                Console.WriteLine("ERROR: SC2 Index out of bounds...");
                return;
            }

            int vectorLength = newColumn.Length;
            if (rows == 0)
            {
                rows = vectorLength;
            }

            if (vectorLength > rows)
            {
                Console.WriteLine("Warning: Input column has too many elements; truncating vector...");
            }
            else if (vectorLength < rows)
            {
                Console.WriteLine("Warning: Input column has too few elements; missing elements will be substituted with 0...");
                double[] subCol = new double[rows];
                for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
                {
                    if (rowIndex < vectorLength)
                    {
                        subCol[rowIndex] = newColumn[rowIndex];
                    }
                    else
                    {
                        subCol[rowIndex] = 0;
                    }
                }
                newColumn = subCol;
            }

            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                values[rowIndex, colIndex] = newColumn[rowIndex];
            }
        }

        public void AddRow(int index, double[] newRow)
        {
            int vectorLength = newRow.Length;
            if (columns == 0)
            {
                columns = vectorLength;
            }

            if (vectorLength > columns)
            {
                Console.WriteLine("Warning: Input row has too many elements; truncating vector...");
            }
            else if (vectorLength < columns)
            {
                Console.WriteLine("Warning: Input row has too few elements; missing elements will be substituted with 0...");
                double[] subRow = new double[columns];
                for (uint colIndex = 0; colIndex < columns; colIndex++)
                {
                    if (colIndex < vectorLength)
                    {
                        subRow[colIndex] = newRow[colIndex];
                    }
                    else
                    {
                        subRow[colIndex] = 0;
                    }
                }
                newRow = subRow;
            }

            if (index > rows)
            {
                Console.WriteLine("Warning: Index input out of bounds; Appending vector to last row of matrix...");
                index = rows;
            }

            rows++;
            double[,] updatedValues = new double[rows, columns];

            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                if (rowIndex < index)
                {
                    for (uint colIndex = 0; colIndex < columns; colIndex++)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex, colIndex];
                    }
                }
                else if (rowIndex > index)
                {
                    for (uint colIndex = 0; colIndex < columns; colIndex++)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex-1, colIndex];
                    }
                }
                else
                {
                    for (uint colIndex = 0; colIndex < columns; colIndex++)
                    {
                        updatedValues[rowIndex, colIndex] = newRow[colIndex];
                    }
                }
            }

            values = updatedValues;
        }
        
        public void RemoveRow(uint index)
        {
            if (index >= rows)
            {
                Console.WriteLine("ERROR: RR Index out of bounds...");
                return;
            }
            rows--;
            double[,] updatedValues = new double[rows, columns];
            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                if (rowIndex < index)
                {
                    for (uint colIndex = 0; colIndex < columns; colIndex++)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex, colIndex];
                    }
                }
                else
                {
                    for (uint colIndex = 0; colIndex < columns; colIndex++)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex+1, colIndex];
                    }
                }
            }
            values = updatedValues;
        }

        public void AddColumn(int index, double[] newColumn)
        {
            int vectorLength = newColumn.Length;
            if (rows == 0)
            {
                rows = vectorLength;
            }

            if (vectorLength > rows)
            {
                Console.WriteLine("Warning: Input column has too many elements; truncating vector...");
            }
            else if (vectorLength < rows)
            {
                Console.WriteLine("Warning: Input column has too few elements; missing elements will be substituted with 0...");
                double[] subCol = new double[rows];
                for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
                {
                    if (rowIndex < vectorLength)
                    {
                        subCol[rowIndex] = newColumn[rowIndex];
                    }
                    else
                    {
                        subCol[rowIndex] = 0;
                    }
                }
                newColumn = subCol;
            }

            if (index > columns)
            {
                Console.WriteLine("Warning: Index input out of bounds; Appending vector to last column of matrix...");
                index = columns;
            }

            columns++;
            double[,] updatedValues = new double[rows, columns];

            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < columns; colIndex++)
                {
                    if (colIndex < index)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex, colIndex];
                    }
                    else if (colIndex > index)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex, colIndex - 1];
                    }
                    else
                    {
                        updatedValues[rowIndex, colIndex] = newColumn[rowIndex];
                    }
                }
            }

            values = updatedValues;
        }

        public void RemoveColumn(uint index)
        {
            if (index >= columns)
            {
                Console.WriteLine("ERROR: RC Index out of bounds...");
                return;
            }
            columns--;
            double[,] updatedValues = new double[rows, columns];
            for (uint rowIndex = 0; rowIndex < rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < columns; colIndex++)
                {
                    if (colIndex < index)
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex, colIndex];
                    }
                    else
                    {
                        updatedValues[rowIndex, colIndex] = values[rowIndex, colIndex+1];
                    }
                }
            }
            values = updatedValues;
        }

        public static double[,] Convert1DArrTo2D(uint inputRows, uint inputCols, double[] inputArray)
        {
            double[,] outputArray = new double[inputRows,inputCols];
            uint arrayIndex = 0;
            int length = inputArray.Length;
            
            if (length < inputRows * inputCols)
            {
                Console.Write("Warning: Not enough values to satisfy ");
                Console.Write(inputRows);
                Console.Write(" by ");
                Console.Write(inputCols);
                Console.WriteLine(" array. Missing values will be substituted with 0... ");
            }
            if (length > inputRows * inputCols)
            {
                Console.Write("Warning: Too many values in array. Array will be truncated to satisfy ");
                Console.Write(inputRows);
                Console.Write(" by ");
                Console.Write(inputCols);
                Console.WriteLine(" array...");
            }

            for (uint rowIndex = 0; rowIndex < inputRows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < inputCols; colIndex++)
                {
                    if (arrayIndex < length)
                    {
                        outputArray[rowIndex, colIndex] = inputArray[arrayIndex];
                    }
                    else
                    {
                        outputArray[rowIndex, colIndex] = 0;
                    }
                    arrayIndex++;
                }
            }
            return outputArray;
        }

        public static double[] AddArrays (double[] array1, double[] array2)
        {
            if (array1.Length > array2.Length)
            {
                double[] sumArray = array1;
                for (uint vectorIndex = 0; vectorIndex < array2.Length; vectorIndex++)
                {
                    sumArray[vectorIndex] += array2[vectorIndex]; 
                }
                return sumArray;
            }
            else
            {
                double[] sumArray = array2;
                for (uint vectorIndex = 0; vectorIndex < array1.Length; vectorIndex++)
                {
                    sumArray[vectorIndex] += array1[vectorIndex];
                }
                return sumArray;
            }
        }

        public static double[] MultiplyArray (double[] array, double scalar)
        {
            double[] product = array;
            for (uint vectorIndex = 0; vectorIndex < array.Length; vectorIndex++)
            {
                product[vectorIndex] = product[vectorIndex] * scalar;
            }
            return product;
        }

        public static double[] DivideArray(double[] array, double scalar)
        {
            double[] quotient = array;
            for (uint vectorIndex = 0; vectorIndex < array.Length; vectorIndex++)
            {
                quotient[vectorIndex] = quotient[vectorIndex] / scalar;
            }
            return quotient;
        }

        public static void Display (double[] vector)
        {
            Console.Write("[ ");
            int vectorLength = vector.Length;
            for (uint vectorIndex = 0; vectorIndex < vectorLength; vectorIndex++)
            {
                Console.Write(vector[vectorIndex] + "\t");
            }
            Console.WriteLine("]");
        }

        public static void Display (Matrix matrix)
        {
            if (matrix == null)
            {
                Console.WriteLine("ERROR: Unable to print matrix...");
                return;
            }

            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                if (rowIndex == 0)
                {
                    Console.Write("┌ ");
                }
                else if (rowIndex == matrix.rows - 1)
                {
                    Console.Write("└ ");
                }
                else
                {
                    Console.Write("| ");
                }

                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    Console.Write(matrix.values[rowIndex, colIndex].ToString("0.###") + "\t");
                }
                if (rowIndex == 0)
                {
                    Console.WriteLine("┐");
                }
                else if (rowIndex == matrix.rows - 1)
                {
                    Console.WriteLine("┘");
                }
                else
                {
                    Console.WriteLine("| ");
                }

            }
        }

        public static Matrix Identity (int dimensions)
        {
            if (dimensions <= 0)
            {
                Console.WriteLine("ERROR: Matrix dimensions must be greater than 0...");
                return null;
            }

            double[,] identityValues = new double[dimensions,dimensions];
            for (uint diagonalIndex = 0; diagonalIndex < dimensions; diagonalIndex++)
            {
                identityValues[diagonalIndex, diagonalIndex] = 1;
            }
            Matrix identity = new Matrix(identityValues);
            return identity;
        }

        public static Matrix I (int dimensions)
        {
            return Identity(dimensions);
        }

        public static Matrix Transpose (Matrix matrix)
        {
            double[,] transposeValues = new double[matrix.columns, matrix.rows];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    transposeValues[colIndex, rowIndex] = matrix.values[rowIndex,colIndex];
                }
            }
            Matrix transpose = new Matrix(transposeValues);
            return transpose;
        }

        public static Matrix T (Matrix matrix)
        {
            return Transpose(matrix);
        }

        public static Matrix operator + (Matrix lhsMatrix, Matrix rhsMatrix)
        {
            if (lhsMatrix.rows != rhsMatrix.rows || lhsMatrix.columns != rhsMatrix.columns)
            {
                Console.WriteLine("ERROR: Cannot sum 2 matrices with different dimensions...");
                return null;
            }

            int sumRows = lhsMatrix.rows;
            int sumCols = lhsMatrix.columns;

            double[,] sumValues = new double[sumRows, sumCols];
            for (uint rowIndex = 0; rowIndex < sumRows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < sumCols; colIndex++)
                {
                    sumValues[rowIndex, colIndex] = lhsMatrix.values[rowIndex, colIndex] + rhsMatrix.values[rowIndex,colIndex];
                }
            }

            Matrix sum = new Matrix(sumValues);
            return sum;
        }

        public static Matrix operator + (Matrix matrix, double scalar)
        {
            double[,] sumValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    sumValues[rowIndex, colIndex] = matrix.values[rowIndex, colIndex] + scalar;
                }
            }
            Matrix sum = new Matrix(sumValues);
            return sum;
        }

        public static Matrix operator + (double scalar, Matrix matrix)
        {
            Matrix sum = matrix + scalar;
            return sum;
        }

        public static Matrix operator - (Matrix matrix)
        {
            double[,] negValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    negValues[rowIndex, colIndex] = -matrix.values[rowIndex, colIndex];
                }
            }
            Matrix negative = new Matrix(negValues);
            return negative;
        }

        public static Matrix operator - (Matrix lhsMatrix, Matrix rhsMatrix)
        {
            if (lhsMatrix.rows != rhsMatrix.rows || lhsMatrix.columns != rhsMatrix.columns)
            {
                Console.WriteLine("ERROR: Cannot sum 2 matrices with different dimensions...");
                return null;
            }

            int difRows = lhsMatrix.rows;
            int difCols = lhsMatrix.columns;

            double[,] difValues = new double[difRows, difCols];
            for (uint rowIndex = 0; rowIndex < difRows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < difCols; colIndex++)
                {
                    difValues[rowIndex, colIndex] = lhsMatrix.values[rowIndex, colIndex] - rhsMatrix.values[rowIndex, colIndex];
                }
            }

            Matrix difference = new Matrix(difValues);
            return difference;
        }

        public static Matrix operator - (Matrix matrix, double scalar)
        {
            double[,] difValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    difValues[rowIndex, colIndex] = matrix.values[rowIndex, colIndex] - scalar;
                }
            }
            Matrix difference = new Matrix(difValues);
            return difference;
        }

        public static Matrix operator - (double scalar, Matrix matrix)
        {
            double[,] difValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    difValues[rowIndex, colIndex] = scalar - matrix.values[rowIndex, colIndex];
                }
            }
            Matrix difference = new Matrix(difValues);
            return difference;
        }

        public static double DotProduct (double[] rowVector, double[] colVector)
        {
            if (rowVector.Length != colVector.Length)
            {
                Console.WriteLine("ERROR: Unequal vector lengths...");
                return 0;
            }
            int vectorLength = rowVector.Length;
            double dotProduct = 0;
            for (uint vectorIndex = 0; vectorIndex < vectorLength; vectorIndex++)
            {
                dotProduct = dotProduct + (rowVector[vectorIndex] * colVector[vectorIndex]);
            }
            return dotProduct;
        }

        public static Matrix operator * (Matrix matrix, double scalar)
        {
            double[,] productValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    productValues[rowIndex, colIndex] = scalar * matrix.values[rowIndex, colIndex];
                }
            }
            Matrix product = new Matrix(productValues);
            return product;
        }

        public static Matrix operator * (double scalar, Matrix matrix)
        {
            Matrix product = matrix * scalar;
            return product;
        }

        public static Matrix operator * (Matrix lhsMatrix, Matrix rhsMatrix)
        {
            if (lhsMatrix.columns != rhsMatrix.rows)
            {
                Console.WriteLine("ERROR: Matrix dimensions incompatible...");
                return null;
            }

            int productRows = lhsMatrix.rows;
            int productCols = rhsMatrix.columns;
            double[,] productValues = new double[productRows, productCols];
            for (uint rowIndex = 0; rowIndex < productRows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < productCols; colIndex++)
                {
                    productValues[rowIndex, colIndex] = DotProduct(lhsMatrix.RowValues(rowIndex), rhsMatrix.ColumnValues(colIndex));
                }
            }
            Matrix product = new Matrix(productValues);
            return product;
        }

        public static Matrix operator ^ (Matrix matrix, uint scalar)
        {
            Matrix product = Copy(matrix);
            for (uint expIndex = 0; expIndex < scalar; expIndex++)
            {
                product = product * matrix;
            }
            return product;
        }

        public static Matrix operator ++ (Matrix matrix)
        {
            double[,] incValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    incValues[rowIndex, colIndex] = matrix.values[rowIndex, colIndex] + 1;
                }
            }
            Matrix incrementedMatrix = new Matrix(incValues);
            return incrementedMatrix;
        }

        public static Matrix operator -- (Matrix matrix)
        {
            double[,] decValues = new double[matrix.rows, matrix.columns];
            for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
            {
                for (uint colIndex = 0; colIndex < matrix.columns; colIndex++)
                {
                    decValues[rowIndex, colIndex] = matrix.values[rowIndex, colIndex] - 1;
                }
            }
            Matrix decrementedMatrix = new Matrix(decValues);
            return decrementedMatrix;
        }

        public static Matrix ReducedEchelonForm (Matrix matrix)
        {
            Matrix ReducedEchelonForm = Copy(matrix);
            uint colIndex = 0;
            uint rowIndex = 0;
            while (colIndex < ReducedEchelonForm.columns && rowIndex < ReducedEchelonForm.rows)
            {
                uint swapIndex = rowIndex + 1;
                while (ReducedEchelonForm.values[rowIndex,colIndex] == 0)
                {
                    if (swapIndex >= ReducedEchelonForm.rows)
                    {
                        return ReducedEchelonForm;
                    }
                    if (ReducedEchelonForm.values[swapIndex,colIndex] != 0)
                    {
                        ReducedEchelonForm.SwapRows(rowIndex, swapIndex);
                    }
                    else
                    {
                        swapIndex++;
                        if (swapIndex >= ReducedEchelonForm.rows)
                        {
                            swapIndex = rowIndex + 1;
                            colIndex++;
                            if (colIndex >= ReducedEchelonForm.columns)
                            {
                                return ReducedEchelonForm;
                            }
                        }
                    }
                }
                double pivot = ReducedEchelonForm.values[rowIndex, colIndex];
                ReducedEchelonForm.SwapRows(rowIndex, DivideArray(ReducedEchelonForm.RowValues(rowIndex), pivot));
                for (uint elimIndex = rowIndex + 1; elimIndex < ReducedEchelonForm.rows; elimIndex++)
                {
                    if (ReducedEchelonForm.values[elimIndex,colIndex] != 0)
                    {
                        ReducedEchelonForm.SwapRows(elimIndex, AddArrays(ReducedEchelonForm.RowValues(elimIndex), MultiplyArray(ReducedEchelonForm.RowValues(rowIndex),-ReducedEchelonForm.values[elimIndex,colIndex])));
                    }
                }
                colIndex++;
                rowIndex++;
            }
            return ReducedEchelonForm;
        }
        public static Matrix REF(Matrix matrix)
        {
            return ReducedEchelonForm(matrix);
        }

        public static Matrix RowReducedEchelonForm(Matrix matrix)
        {
            Matrix RowReducedEchelonForm = Copy(matrix);
            uint colIndex = 0;
            uint rowIndex = 0;
            while (colIndex < RowReducedEchelonForm.columns && rowIndex < RowReducedEchelonForm.rows)
            {
                uint swapIndex = rowIndex + 1;
                while (RowReducedEchelonForm.values[rowIndex, colIndex] == 0)
                {
                    if (swapIndex >= RowReducedEchelonForm.rows)
                    {
                        return RowReducedEchelonForm;
                    }
                    if (RowReducedEchelonForm.values[swapIndex, colIndex] != 0)
                    {
                        RowReducedEchelonForm.SwapRows(rowIndex, swapIndex);
                    }
                    else
                    {
                        swapIndex++;
                        if (swapIndex >= RowReducedEchelonForm.rows)
                        {
                            swapIndex = rowIndex + 1;
                            colIndex++;
                            if (colIndex >= RowReducedEchelonForm.columns)
                            {
                                return RowReducedEchelonForm;
                            }
                        }
                    }
                }
                double pivot = RowReducedEchelonForm.values[rowIndex, colIndex];
                RowReducedEchelonForm.SwapRows(rowIndex, DivideArray(RowReducedEchelonForm.RowValues(rowIndex), pivot));
                for (uint elimIndex = 0; elimIndex < RowReducedEchelonForm.rows; elimIndex++)
                {
                    if (RowReducedEchelonForm.values[elimIndex, colIndex] != 0 && elimIndex != rowIndex)
                    {
                        RowReducedEchelonForm.SwapRows(elimIndex, AddArrays(RowReducedEchelonForm.RowValues(elimIndex), MultiplyArray(RowReducedEchelonForm.RowValues(rowIndex), -RowReducedEchelonForm.values[elimIndex, colIndex])));
                    }
                }
                colIndex++;
                rowIndex++;
            }
            return RowReducedEchelonForm;
        }
        public static Matrix RREF (Matrix matrix)
        {
            return RowReducedEchelonForm(matrix);
        }

        public static (Matrix, Matrix, Matrix) PLUDecomposition (Matrix matrix)
        {
            Matrix PermutationMatrix = Identity(matrix.rows);
            Matrix LowerTriangularMatrix;
            if (matrix.rows > matrix.columns)
            {
                LowerTriangularMatrix = new Matrix(matrix.rows, matrix.columns);
            }
            else
            {
                LowerTriangularMatrix = new Matrix(matrix.rows, matrix.rows);
            }
            Matrix UpperTriangularMatrix = Copy(matrix);

            for (uint colIndex = 0; colIndex < LowerTriangularMatrix.columns; colIndex++)
            {
                uint pivotIndex = colIndex;
                double pivotValue = UpperTriangularMatrix.values[colIndex, colIndex];
                if (pivotValue < 0)
                {
                    pivotValue = pivotValue * -1;
                }
                for (uint rowIndex = colIndex + 1; rowIndex < LowerTriangularMatrix.rows; rowIndex++)
                {
                    double rowValue = UpperTriangularMatrix.values[rowIndex, colIndex];
                    if (rowValue < 0)
                    {
                        rowValue = rowValue * -1;
                    }
                    if (pivotValue < rowValue)
                    {
                        pivotIndex = rowIndex;
                        pivotValue = rowValue;
                    }
                }
                if (pivotIndex != colIndex)
                {
                    PermutationMatrix.SwapRows(colIndex, pivotIndex);
                    LowerTriangularMatrix.SwapRows(colIndex, pivotIndex);
                    UpperTriangularMatrix.SwapRows(colIndex, pivotIndex);
                }

                LowerTriangularMatrix.values[colIndex, colIndex] = 1;
                if (UpperTriangularMatrix.values[colIndex,colIndex] == 0)
                {
                    continue;
                }

                for (uint rowIndex = colIndex + 1; rowIndex < LowerTriangularMatrix.rows; rowIndex++)
                {
                    double operation = -(UpperTriangularMatrix.values[rowIndex, colIndex] / UpperTriangularMatrix.values[colIndex, colIndex]);
                    LowerTriangularMatrix.values[rowIndex, colIndex] = -operation;
                    UpperTriangularMatrix.SwapRows(rowIndex, AddArrays(MultiplyArray(UpperTriangularMatrix.RowValues(colIndex), operation), UpperTriangularMatrix.RowValues(rowIndex)));
                }
            }

            for (uint rowIndex = (uint)matrix.columns; rowIndex < matrix.rows; rowIndex++)
            {
                UpperTriangularMatrix.RemoveRow((uint)matrix.columns);
            }
            PermutationMatrix = Transpose(PermutationMatrix);
            return (PermutationMatrix, LowerTriangularMatrix, UpperTriangularMatrix);
        }

        public static (Matrix, Matrix, Matrix) PLU (Matrix matrix)
        {
            return PLUDecomposition(matrix);
        }

        public static Matrix Inverse (Matrix matrix)
        {
            if (matrix.rows != matrix.columns)
            {
                Console.WriteLine("ERROR: Non-square matrix not invertible...");
                return null;
            }

            Matrix augment = Copy(matrix);
            Matrix inverse = Identity(matrix.rows);

            for (uint diagIndex = 0; diagIndex < matrix.rows; diagIndex++)
            {
                if (augment.values[diagIndex,diagIndex] == 0)
                {
                    Console.WriteLine("ERROR: Matrix not invertible...");
                    return null;
                }
                inverse.SwapRows(diagIndex, DivideArray(inverse.RowValues(diagIndex), augment.values[diagIndex, diagIndex]));
                augment.SwapRows(diagIndex, DivideArray(augment.RowValues(diagIndex), augment.values[diagIndex, diagIndex]));
                for (uint rowIndex = 0; rowIndex < matrix.rows; rowIndex++)
                {
                    if(rowIndex == diagIndex)
                    {
                        continue;
                    }
                    inverse.SwapRows(rowIndex, AddArrays(inverse.RowValues(rowIndex), MultiplyArray(inverse.RowValues(diagIndex), -augment.values[rowIndex, diagIndex])));
                    augment.SwapRows(rowIndex, AddArrays(augment.RowValues(rowIndex), (MultiplyArray(augment.RowValues(diagIndex), -augment.values[rowIndex,diagIndex]))));
                }
            }
            return inverse;
        }
    }   
}
