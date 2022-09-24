#include "SelectionParser.hpp"
#include "selections/selections.hpp"
#include <stack>
#include <string_view>
#include <molpp/MolError.hpp>

using namespace mol;
using namespace mol::internal;

struct CopyNode : public SelectionNode
{
    CopyNode(std::shared_ptr<peg::Ast> const& node)
    : ast(node)
    {}

    void evaluate(SelectionStack& /*evaluator*/, MolData const& /*data*/, Frame /*frame*/) const override
    {};

    std::shared_ptr<peg::Ast> ast;
};

std::shared_ptr<SelectionNode> make_boolbinary_node(std::shared_ptr<peg::Ast> const ast)
{
    std::shared_ptr<SelectionNode> node;
    if (ast->nodes[1]->token == "and")
    {
        node = std::make_shared<AndSelection>();
    }
    else if (ast->nodes[1]->token == "or")
    {
        node = std::make_shared<OrSelection>();
    }

    // Set temporary children holding the Ast nodes
    node->left = std::make_shared<CopyNode>(ast->nodes[0]);
    node->right = std::make_shared<CopyNode>(ast->nodes[2]);

    return node;
}

std::shared_ptr<SelectionNode> make_boolunary_node(std::shared_ptr<peg::Ast> const ast)
{
    std::shared_ptr<SelectionNode> node = std::make_shared<NotSelection>();
    // Set temporary children holding the Ast nodes
    node->left = std::make_shared<CopyNode>(ast->nodes[1]);
    return node;
}

std::shared_ptr<SelectionNode> make_numprop_node(std::shared_ptr<peg::Ast> const ast)
{
    std::shared_ptr<NumPropSelection> num_prop;
    auto const& type = ast->nodes[0]->token;
    if (type == "resid")
    {
        num_prop = std::make_shared<ResidSelection>();
    }
    else
    {
        // We should never get here
        throw mol::MolError("Unknown NumProp node: " + std::string(type));
    }
    std::shared_ptr<SelectionNode> node = num_prop;

    for (std::shared_ptr<peg::Ast> child : ast->nodes)
    {
        if (child->name == "Number")
        {
            std::string_view const& token = child->token;
            num_prop->add_number(SelNumber({token.data(), token.size()}));
        }
        else if (child->name == "NumRange")
        {
            std::string_view const& token1 = child->nodes[0]->token;
            std::string_view const& token2 = child->nodes[1]->token;
            num_prop->add_range(SelNumberRange({token1.data(), token1.size()}, {token2.data(), token2.size()}));
        }
    }

    return node;
}

std::shared_ptr<SelectionNode> make_node(std::shared_ptr<peg::Ast> const ast)
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

    return nullptr;
}

SelectionParser::SelectionParser(std::string const& grammar)
: m_grammar{grammar}
{
    if (!m_parser.load_grammar(m_grammar))
    {
        throw mol::MolError("Error loading selection grammar.");
    }

    m_parser.log = [&](size_t /*line*/, size_t column, std::string const& message)
    {
        m_error_column = column;
        m_error_message = message;
    };
    m_parser.enable_ast();
}

std::shared_ptr<SelectionNode> SelectionParser::parse(std::string const& expression) const
{
    // Parsing
    std::shared_ptr<peg::Ast> ast;
    if (!m_parser.parse(expression, ast))
    {
        throw mol::MolError(error_message(expression));
    }
    ast = m_parser.optimize_ast(ast);

    // Convertion into selection tree
    std::stack<std::shared_ptr<SelectionNode>> node_stack;
    std::shared_ptr<SelectionNode> root = make_node(ast);
    if (!root)
    {
        throw mol::MolError("Error while building the selection tree.");
    }

    node_stack.push(root);
    while (!node_stack.empty())
    {
        std::shared_ptr<SelectionNode> current = node_stack.top();
        node_stack.pop();

        // Right child
        if (current->right)
        {
            std::shared_ptr<peg::Ast> right_ast = std::static_pointer_cast<CopyNode>(current->right)->ast;
            std::shared_ptr<SelectionNode> right = make_node(right_ast);
            if (!right)
            {
                throw mol::MolError("Error while building the selection tree.");
            }
            current->right = right;
            node_stack.push(right);
        }

        // Left child
        if (current->left)
        {
            std::shared_ptr<peg::Ast> left_ast = std::static_pointer_cast<CopyNode>(current->left)->ast;
            std::shared_ptr<SelectionNode> left = make_node(left_ast);
            if (!left)
            {
                throw mol::MolError("Error while building the selection tree.");
            }
            current->left = left;
            node_stack.push(left);
        }
    }

    return root;
}

std::string SelectionParser::error_message(std::string const& expression) const
{
    std::string error = m_error_message + "\n" + expression + "\n";
    for (size_t i = 1; i < m_error_column; i++)
    {
        error += " ";
    }
    error += "^";

    return error;
}

extern SelectionParser const mol::internal::SEL_PARSER(R"(
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
