#include <iostream>
using namespace std;
class test
{
private:
	double a[4] = {1.0,2.0,3.0,4.0};
public:
	test()
	{
		//*a = {1,2,3,4};
	}
	void pf()
	{
		for(int i = 0; i < 4; ++i)
		{
			cout<<a[i]<<endl;
		}
	}
};

int main()
{
	test * p = new test;
	p->pf();
	return 0;
}

