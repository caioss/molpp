#ifndef MOLPP_SELECTIONS_SELECTIONPARSER_HPP
#define MOLPP_SELECTIONS_SELECTIONPARSER_HPP

#include <memory>
#include <string>

#include <peglib.h>

namespace mol::internal
{

class SelectionNode;

class SelectionParser
{
public:
    SelectionParser(std::string const& grammar);
    std::shared_ptr<SelectionNode> parse(std::string const& expression);

private:
    std::string m_grammar;
    peg::parser m_parser;
};

SelectionParser& default_parser();

} // namespace mol::internal

#endif // MOLPP_SELECTIONS_SELECTIONPARSER_HPP
