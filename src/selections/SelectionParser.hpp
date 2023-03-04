#ifndef SELECTIONPARSER_HPP
#define SELECTIONPARSER_HPP

#include <memory>
#include <string>
#include <peglib.h>

namespace mol::internal {

class SelectionNode;

class SelectionParser
{
public:
    SelectionParser(std::string const& grammar);
    std::shared_ptr<SelectionNode> parse(std::string const& expression) const;
    std::string error_message(std::string const& expression) const;

private:
    std::size_t m_error_column;
    std::string m_grammar;
    std::string m_error_message;
    peg::parser m_parser;
};

extern SelectionParser const SEL_PARSER;

} // namespace mol::internal

#endif // SELECTIONPARSER_HPP
