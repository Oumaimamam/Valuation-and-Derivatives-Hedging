#include <iostream>
#include <fstream>
#include "../Utils/convert.cpp"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Taille attendu 2" << endl;
        return 1;
    }

    fstream file(argv[1]);
    nlohmann::json json = nlohmann::json::parse(file);

    Option option;

    convert_json_to_option(&option, json);

    return 0;
}
