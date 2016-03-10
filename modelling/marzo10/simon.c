#include <stdio.h>

void fill(int n, int *v)
{
  int i,x;
  for (i=0; i<n; i=i+1)
{
    scanf("%d", &x);
    v[i]= x;
 }
}

int main()
{
  int i, v[5], x;
  fill(5, &v);
  for (i=0; i<5;i=i+1) printf("%d ",v[i]);
  system("PAUSE");
  return (0);
}
