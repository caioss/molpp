#include "SelectionParser.hpp"
#include "selections/boolean.hpp"
#include "selections/properties.hpp"
#include "selections/NumberSet.hpp"
#include <molpp/Error.hpp>

#include <stack>
#include <ranges>
#include <string_view>

namespace mol::internal
{

class SelectionStack;

//! Helper node that holds a copy of the AST node for later use.
struct CopyNode : public SelectionNode
{
    CopyNode(std::shared_ptr<peg::Ast> const& node)
    : ast(node)
    {}

    void evaluate(SelectionStack& /*evaluator*/, MolData const& /*data*/, Frame /*frame*/) const override
    {};

    std::shared_ptr<peg::Ast> ast;
};

std::shared_ptr<SelectionNode> make_boolbinary_node(std::shared_ptr<peg::Ast const> const ast)
{
    std::shared_ptr<SelectionNode> node;
    std::string_view const type = ast->nodes[1]->token;
    if (type == "and")
    {
        node = std::make_shared<AndSelection>();
    }
    else if (type == "or")
    {
        node = std::make_shared<OrSelection>();
    }
    else
    {
        // We should never get here
        throw mol::Error("Unknown boolean node: " + std::string(type));
    }

    // Set temporary children holding the AST nodes
    node->left = std::make_shared<CopyNode>(ast->nodes[0]);
    node->right = std::make_shared<CopyNode>(ast->nodes[2]);

    return node;
}

std::shared_ptr<SelectionNode> make_boolunary_node(std::shared_ptr<peg::Ast const> const ast)
{
    std::shared_ptr<SelectionNode> node = std::make_shared<NotSelection>();
    // Set temporary children holding the AST nodes
    node->left = std::make_shared<CopyNode>(ast->nodes[1]);
    return node;
}

template<class Type>
std::shared_ptr<SelectionNode> make_numprop_node_impl(std::vector<std::shared_ptr<peg::Ast>> const& numbers)
{
    NumberSet number_set;

    // Skip the first node, which is the property node
    for (std::shared_ptr<const peg::Ast> child : numbers | std::views::drop(1))
    {
        std::string const& type = child->name;
        if (type == "Number")
        {
            std::string const token{child->token};
            number_set.add_number(SelNumber(token));
        }
        else if (type == "NumRange")
        {
            std::string const token1{child->nodes[0]->token};
            std::string const token2{child->nodes[1]->token};
            number_set.add_range(SelNumberRange(token1, token2));
        }
        else
        {
            // We should never get here
            throw mol::Error("Unknown numeric node: " + std::string(type));
        }
    }

    std::shared_ptr<SelectionNode> node = std::make_shared<Type>(std::move(number_set));
    return node;
}

std::shared_ptr<SelectionNode> make_numprop_node(std::shared_ptr<peg::Ast const> const ast)
{
    std::string_view const type = ast->nodes[0]->token;
    if (type == "resid")
    {
        return make_numprop_node_impl<ResidSelection>(ast->nodes);
    }
    else
    {
        // We should never get here
        throw mol::Error("Unknown NumProp node: " + std::string(type));
    }
}

std::shared_ptr<SelectionNode> make_all_node(std::shared_ptr<peg::Ast const> const /*ast*/)
{
    return std::make_shared<AllSelection>();
}

std::shared_ptr<SelectionNode> make_node(std::shared_ptr<peg::Ast const> const ast)
{
    /* Boolean binary operators */
    if (ast->name == "BoolBinaryExp")
    {
        return make_boolbinary_node(ast);
    }

    /* Boolean unary */
    else if (ast->name == "BoolUnaryExp")
    {
        return make_boolunary_node(ast);
    }

    /* Numerical atom properties */
    else if (ast->name == "NumPropExp")
    {
        return make_numprop_node(ast);
    }

    /* All */
    else if (ast->name == "All")
    {
        return make_all_node(ast);
    }

    return nullptr;
}

std::shared_ptr<SelectionNode> make_child(std::shared_ptr<SelectionNode> const node)
{
    std::shared_ptr<peg::Ast> ast = std::static_pointer_cast<CopyNode>(node)->ast;
    std::shared_ptr<SelectionNode> child = make_node(ast);
    if (!child)
    {
        throw mol::Error("Error while building the selection tree.");
    }

    return child;
}

class ParserError
{
public:
    ParserError(std::string const& expression)
    : m_expression{expression}
    {}

    void set_error(size_t column, std::string const& message)
    {
        m_error_column = column;
        m_error_message = message;
    }

    std::string message() const
    {
        std::string error = m_error_message + "\n" + m_expression + "\n";
        for (size_t i = 1; i < m_error_column; i++)
        {
            error += " ";
        }
        error += "^";

        return error;
    }

private:
    std::size_t m_error_column;
    std::string m_error_message;
    std::string const& m_expression;
};

SelectionParser::SelectionParser(std::string const& grammar)
: m_grammar{grammar}
{
    if (!m_parser.load_grammar(m_grammar))
    {
        throw mol::Error("Error loading grammar.");
    }

    m_parser.enable_ast();
}

std::shared_ptr<SelectionNode> SelectionParser::parse(std::string const& expression)
{
    // Parsing
    ParserError error(expression);
    m_parser.set_logger([&](size_t /*line*/, size_t column, std::string const& message) {
        error.set_error(column, message);
    });

    std::shared_ptr<peg::Ast> ast;
    if (!m_parser.parse(expression, ast))
    {
        throw mol::Error(error.message());
    }
    ast = m_parser.optimize_ast(ast);

    // Convert AST into a selection tree.
    // We don't use recursion to allow huge expressions.
    std::stack<std::shared_ptr<SelectionNode>> node_stack;
    std::shared_ptr<SelectionNode> root = make_node(ast);
    if (!root)
    {
        throw mol::Error("Error while building the selection tree.");
    }

    node_stack.push(root);
    while (!node_stack.empty())
    {
        std::shared_ptr<SelectionNode> current = node_stack.top();
        node_stack.pop();

        if (current->right)
        {
            current->right = make_child(current->right);
            node_stack.push(current->right);
        }

        if (current->left)
        {
            current->left = make_child(current->left);
            node_stack.push(current->left);
        }
    }

    return root;
}

SelectionParser& default_parser()
{
    static SelectionParser parser(R"(
    BoolBinaryExp <- Operand (BoolBinaryOp Operand)* {
                       precedence
                         L and or
                     }
    BoolUnaryExp  <- BoolUnaryOp Operand
    Operand       <- BoolUnaryExp / '(' BoolBinaryExp ')' / NumPropExp / All

    NumPropExp    <- NumProp (NumRange / Number)+
    TextPropExp   <- TextProp Text+
    SameExp       <- 'same' Prop 'as' Operand
    WithinExp     <- 'within' Number 'of' (Operand / Vector)
    CenterExp     <- 'center of' Operand

    All           <- 'all'
    BoolOperators <- BoolUnaryOp / BoolBinaryOp / All
    BoolBinaryOp  <- 'and' / 'or'
    BoolUnaryOp   <- 'not'

    CompOperator  <- '==' / '!=' / '<' / '<=' / '>' / '>='

    NumComp       <- NumExp (CompOperator NumExp){1,2}
    NumExp        <- NumOperand (NumOperator NumOperand)* {
                       precedence
                         L + -
                         L * /
                         R **
                     }
    NumOperator   <-  '-' / '+' / '/' / '**' / '*'
    NumOperand    <-  Number / NumProp / '(' NumExp ')'

    Prop          <- NumProp / TextProp
    NumProp       <- 'resid' / 'index'
    TextProp      <- 'resname' / 'name'

    NumRange      <-  Number ':' Number
    Number        <- < ('-' / '+')? [0-9]+ ('.' [0-9]+ )? >
    Text          <- < !BoolOperators [a-zA-Z0-9_'"]+ >
    Vector        <- '(' Number{3} ')' / CenterExp

    %whitespace   <- [ \t\n]*
    )");

    return parser;
}

} // namespace mol::internal
