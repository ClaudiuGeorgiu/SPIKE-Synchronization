#include "SPIKESynchronization.h"

int main(void)
{
	SPIKESynchronization* spike = new SPIKESynchronization();

	if (spike != NULL)
		delete spike;

	return 0;
}