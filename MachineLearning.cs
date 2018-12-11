using System;
using System.Collections.Generic;
using System.Text;

// MeanSquaredError
// StepGradient
// LinearRegressionGradientDescent (Univariate)
// LinearRegressionNormalEquation (Multivariate)

namespace Matrix
{
    class MachineLearning
    {
        // MeanSquaredError : Returns the mean squared error given a slope, an offset, and a set of points represented by a matrix where column 0 are x coordinates and column 1 are y coordinates
        public static double MeanSquaredError (double offset, double slope, Matrix points)
        {
            double[,] XYpoints = points.Values();
            double error = 0;
            if (points.Columns() != 2)
            {
                Console.WriteLine("ERROR: Mean Squared Error can only handle 2d points...");
                return 0;
            }
            for (uint rowIndex = 0; rowIndex < points.Rows(); rowIndex++)
            {
                error += ((XYpoints[rowIndex, 1] - (slope * XYpoints[rowIndex,0] + offset)) * (XYpoints[rowIndex, 1] - (slope * XYpoints[rowIndex, 0] + offset)));
            }
            return (error/points.Rows());
        }

        // StepGradient : Computes the gradient through the differentiation of the Mean Squared Error function
        public static (double,double) StepGradient(Matrix points, double learningRate, double currentOffset, double currentSlope)
        {
            if (points.Columns() != 2)
            {
                Console.WriteLine("ERROR: Mean Squared Error can only handle 2d points...");
                return (0,0);
            }
            double offsetGradient = 0;
            double slopeGradient = 0;
            double[,] XYpoints = points.Values();
            double numberOfPoints = points.Rows();
            for (uint rowIndex = 0; rowIndex < numberOfPoints; rowIndex++)
            {
                offsetGradient += -(2 / numberOfPoints) * (XYpoints[rowIndex,1] - ((currentSlope * XYpoints[rowIndex,0]) + currentOffset));
                slopeGradient += -(2 / numberOfPoints) * XYpoints[rowIndex, 0] * (XYpoints[rowIndex,1]- ((currentSlope * XYpoints[rowIndex, 0]) + currentOffset));
            }
            double newOffset = currentOffset - (learningRate * offsetGradient);
            double newSlope = currentSlope - (learningRate * slopeGradient);
            return (newOffset, newSlope);
        }

        // StepGradient (Multivariate) : Computes the gradient of a multivariate design matrix and an output
        public static double[] StepGradient (Matrix design, Matrix output, double[] currentWeights, double learningRate)
        {
            double numberOfSamples = design.Rows();
            double numberOfFeatures = design.Columns();
            double[,] designValues = design.Values();
            double[,] outputValues = output.Values();
            double[] gradients = new double[currentWeights.Length];
            if (output.Columns()!=1)
            {
                Console.WriteLine("ERROR: Output vector must be single dimensional...");
                return null;
            }
            if (output.Rows() != numberOfSamples)
            {
                Console.WriteLine("ERROR: Output vector does not have same dimensions as training set...");
                return null;
            }
            if (currentWeights.Length != numberOfFeatures)
            {
                Console.WriteLine("ERROR: Weight vector does not have same length as features...");
                return null;
            }
            for (uint rowIndex = 0; rowIndex < numberOfSamples; rowIndex++)
            {
                double hypothesis = 0;
                for (uint weightIndex = 0; weightIndex < numberOfFeatures; weightIndex++)
                {
                    hypothesis += currentWeights[weightIndex] * designValues[rowIndex, weightIndex];
                }
                for (uint featureIndex = 0; featureIndex < numberOfFeatures; featureIndex++)
                {
                    gradients[featureIndex] += -(2 / numberOfSamples) * designValues[rowIndex, featureIndex] * (outputValues[rowIndex, 0] - hypothesis);
                }
            }
            double[] newWeights = new double[currentWeights.Length];
            for (uint weightIndex = 0; weightIndex < numberOfFeatures; weightIndex++)
            {
                newWeights[weightIndex] = currentWeights[weightIndex] - (learningRate * gradients[weightIndex]);
            }
            return newWeights;
        }

        // LinearRegressionGradientDescent : Finds the closest line that will fit a given set of points from a Matrix (column[0] is X-coordinate and column[1] is Y-coordinate) through gradient descent
        public static Matrix LinearRegressionGradientDescent (Matrix points, double learningRate, double initialOffset, double initialSlope, int iterationLimit, double convergenceLimit)
        {
            if (points.Columns() != 2)
            {
                Console.WriteLine("ERROR: Mean Squared Error can only handle 2d points...");
                return null;
            }
            double offset = initialOffset;
            double slope = initialSlope;
            double previousOffsetShift = double.MaxValue;
            double previousSlopeShift = double.MaxValue;
            int iterations = 0;
            while (iterations < iterationLimit)
            {
                (double,double) gradient = StepGradient(points, learningRate, offset, slope);
                double offsetShift = gradient.Item1 - offset;
                if (offsetShift < 0)
                {
                    offsetShift *= -1;
                }
                double slopeShift = gradient.Item2 - slope;
                if (slopeShift < 0)
                {
                    slopeShift *= -1;
                }
                if ((offsetShift <= convergenceLimit) && (slopeShift <= convergenceLimit))
                {
                    break;
                }
                iterations++;
                if (offsetShift > previousOffsetShift && slopeShift > previousSlopeShift)
                {
                    learningRate *= 0.1;
                    continue;
                }
                offset = gradient.Item1;
                slope = gradient.Item2;
                previousOffsetShift = offsetShift;
                previousSlopeShift = slopeShift;
            }
            Console.WriteLine("Gradient Descent completed in {0} iterations...", iterations);
            Console.WriteLine("Learning Rate: {0}", learningRate);
            Console.WriteLine("Mean Squared Error = {0}", MeanSquaredError(offset, slope, points));
            double[,] vector = { { offset },{ slope } };
            Matrix solution = new Matrix(vector);
            return solution;
        }

        // LinearRegressionGradient (Multivariate) : Calculates the theta vector to determine the linear regression of a multivariate system 
        // TO DO: Optimize
        public static Matrix LinearRegressionGradientDescent (Matrix design, Matrix output, double[] initialWeights, double learningRate, int iterationLimit, double convergenceLimit)
        {
            if (output.Columns() != 1)
            {
                Console.WriteLine("ERROR: Output vector must be single dimensional...");
                return null;
            }
            if (output.Rows() != design.Rows())
            {
                Console.WriteLine("ERROR: Output vector does not have same dimensions as training set...");
                return null;
            }

            Matrix input = Matrix.Copy(design);
            double[] theta0 = new double[design.Rows()];
            for (uint rowIndex = 0; rowIndex < design.Rows(); rowIndex++)
            {
                theta0[rowIndex] = 1;
            }
            input.AddColumn(0, theta0);
            if (initialWeights.Length != input.Columns())
            {
                Console.WriteLine("ERROR: Weight vector does not have same length as features...");
                return null;
            }

            int numberOfFeatures = input.Columns();
            double[] weights = new double[numberOfFeatures];
            for (uint weightIndex = 0; weightIndex < numberOfFeatures; weightIndex++)
            {
                weights[weightIndex] = initialWeights[weightIndex];
            }
            double prevShift = double.MaxValue;
            uint iterations = 0;
            while (iterations < iterationLimit)
            {
                iterations++;
                double[] newWeights = StepGradient(input, output, weights, learningRate);
                double[] shift = new double[numberOfFeatures];
                int convergenceCheck = 0;
                for (uint shiftIndex = 0; shiftIndex < numberOfFeatures; shiftIndex++)
                {
                    double offsetShift = weights[shiftIndex] - newWeights[shiftIndex];
                    if (offsetShift < 0)
                    {
                        offsetShift *= -1;
                    }
                    if (offsetShift <= convergenceLimit)
                    {
                        convergenceCheck++;
                    }
                    shift[shiftIndex] = offsetShift;
                }
                if (convergenceCheck == numberOfFeatures)
                {
                    break;
                }
                double currentShift = Matrix.VectorLength(shift);
                if (currentShift > prevShift)
                {
                    learningRate = (learningRate * 0.1);
                    continue;
                }
                prevShift = currentShift;
                for (uint weightIndex = 0; weightIndex < numberOfFeatures; weightIndex++)
                {
                    weights[weightIndex] = newWeights[weightIndex];
                }
            }
            double[,] vector = Matrix.Convert1DArrTo2D((uint)numberOfFeatures, 1, weights);
            Console.WriteLine("Gradient Descent completed in {0} iterations...", iterations);
            Console.WriteLine("Learning Rate: {0}", learningRate);
            Matrix solution = new Matrix(vector);
            return solution;
        }

        // LinearRegressionNormalEquation: Finds the closest line that will fit a given design matrix with multiple features and an output vector as a result of the design matrix
        public static Matrix LinearRegressionNormalEquation (Matrix design, Matrix output)
        {
            if (output.Columns() != 1)
            {
                Console.WriteLine("ERROR: Output vector not unidimensional...");
                return null;
            }
            Matrix input = Matrix.Copy(design);
            double[] theta0 = new double[design.Rows()];
            for (uint rowIndex = 0; rowIndex < design.Rows(); rowIndex++)
            {
                theta0[rowIndex] = 1;
            }
            input.AddColumn(0, theta0);
            Matrix weight = ((Matrix.T(input)*input)^-1)*(Matrix.T(input)*output);
            return weight;
        }
    }
}
