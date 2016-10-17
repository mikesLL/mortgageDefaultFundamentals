# mortgageDefaultFundamentals
This project models mortgage default in a household's life-cycle portfolio optimization problem. The model incorporates a home price process which displays short-run serial correlation and a long-run relationship with rents. In each period, the household chooses selects the optimal mortgage size, type, and whether to default. In the event of a moving shock, income shock, or decline in home prices, the household may also consider the landlord option, or renting out the house instead of defaulting. 


## Code

main.cpp is the main file in this project.

## Compile Instructions
On Linux:
g++ *.cpp -std=c++11 -fopenmp -lpthread

