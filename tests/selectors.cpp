#include <molpp/AtomSelector.hpp>
#include "selections/SelectionParser.hpp"
#include "selections/selections.hpp"
#include "auxiliary.hpp"
#include "files.hpp"
#include <molpp/MolError.hpp>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <set>

using namespace mol;
using namespace mol::internal;
using namespace testing;

std::shared_ptr<std::set<size_t>> evaluate_sel_tree(std::shared_ptr<SelectionNode> tree, MolData const& data)
{
    SelectionFlags flags;
    for (size_t atom_idx = 0; atom_idx < data.size(); atom_idx++)
    {
        flags.mask->insert(atom_idx);
    }

    SelectionStack eval(tree);
    eval.evaluate(data, flags, 0);
    return flags.selected;
}

TEST(Selection, DataNodes) {
    SelNumber float_num("1.0001");
    EXPECT_EQ(float_num, 1.0001);
    EXPECT_EQ(float_num, 1.0002);
    EXPECT_NE(float_num, 1.0003);

    SelNumber int_num("1");
    EXPECT_EQ(int_num, 1);
    EXPECT_NE(int_num, 2);

    SelNumberRange range("3", "5");
    // Integer
    EXPECT_FALSE(range.has(-1));
    EXPECT_FALSE(range.has(2));
    EXPECT_TRUE(range.has(3));
    EXPECT_TRUE(range.has(4));
    EXPECT_TRUE(range.has(5));
    EXPECT_FALSE(range.has(6));
    // Floating point
    EXPECT_FALSE(range.has(-1.0));
    EXPECT_FALSE(range.has(2.0));
    EXPECT_TRUE(range.has(3.0));
    EXPECT_TRUE(range.has(4.0));
    EXPECT_TRUE(range.has(5.0));
    EXPECT_FALSE(range.has(6.0));

    SelNumberRange consecutive_range("3", "4");
    EXPECT_FALSE(consecutive_range.has(2));
    EXPECT_TRUE(consecutive_range.has(3));
    EXPECT_TRUE(consecutive_range.has(4));
    EXPECT_FALSE(consecutive_range.has(5));
}

TEST(Selection, PropNodes) {
    // NumPropSelection (base)
    class SubNumProp : public NumPropSelection
    {
    public:
        void evaluate(SelectionStack& evaluator, MolData const& data, std::optional<size_t> frame) const override {};
    };

    SubNumProp num_prop;
    EXPECT_FALSE(num_prop.has(1));
    EXPECT_FALSE(num_prop.has(1.0));

    num_prop.add_number(SelNumber("0"));
    num_prop.add_number(SelNumber("2"));
    num_prop.add_range(SelNumberRange("4", "6"));
    num_prop.add_range(SelNumberRange("8", "10"));
    for (int value : {0, 2, 4, 5, 6, 8, 9, 10})
    {
        EXPECT_TRUE(num_prop.has(value)) << value;
        EXPECT_TRUE(num_prop.has(double(value))) << value;
    }
    for (int value : {1, 3, 7, 11})
    {
        EXPECT_FALSE(num_prop.has(value)) << value;
        EXPECT_FALSE(num_prop.has(double(value))) << value;
    }
}

TEST(Selection, ResidSelection) {
    auto data = create_moldata(5, 2, 2, 3, 1);
    SelectionStack flag_stack(nullptr);
    SelectionFlags flags;
    for (size_t atom_idx = 0; atom_idx < data->size(); atom_idx++)
    {
        flags.mask->insert(atom_idx);
    }

    // Default resid
    ResidSelection resid;
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre());

    // One number
    resid.add_number(SelNumber("0"));
    flags.selected->clear();
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre(0, 1));

    // Two numbers
    resid.add_number(SelNumber("4"));
    flags.selected->clear();
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre(0, 1, 8, 9));

    // De-select one atom
    flags.mask->erase(0);
    flags.selected->clear();
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre(1, 8, 9));
    // Bring back mask state
    flags.mask->insert(0);

    // One range
    resid.add_range(SelNumberRange("1", "2"));
    flags.selected->clear();
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre(0, 1, 2, 3, 4, 5, 8, 9));

    // Two ranges
    resid.add_range(SelNumberRange("1", "3"));
    flags.selected->clear();
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre(0, 1, 2, 3, 4, 5, 6, 7, 8, 9));

    // De-select some atoms
    flags.mask->erase(0);
    flags.mask->erase(2);
    flags.mask->erase(5);
    flags.mask->erase(7);
    flags.mask->erase(8);
    flags.selected->clear();
    flag_stack.push_flags(flags);
    resid.evaluate(flag_stack, *data, 0);
    EXPECT_THAT(*(flags.selected), UnorderedElementsAre(1, 3, 4, 6, 9));
    // Bring back mask state
    flags.mask->insert(0);
    flags.mask->insert(2);
    flags.mask->insert(5);
    flags.mask->insert(7);
    flags.mask->insert(8);
}

TEST(Selection, Evaluation) {
    // SelectionFlags default constructor
    SelectionFlags flags;
    ASSERT_THAT(flags.mask, NotNull());
    EXPECT_THAT(*(flags.mask), ElementsAre());
    ASSERT_THAT(flags.selected, NotNull());
    EXPECT_THAT(*(flags.selected), ElementsAre());

    // Check if nodes get eventually evaluated
    class DummyNode : public SelectionNode
    {
    public:
        MOCK_METHOD(void, evaluate, (SelectionStack& evaluator, MolData const& data, std::optional<size_t> frame), (const, override));
    };

    std::shared_ptr<SelectionNode> root = std::make_shared<AndSelection>();
    auto left = std::make_shared<DummyNode>();
    root->left = left;
    EXPECT_CALL(*left, evaluate).Times(1);
    auto right = std::make_shared<DummyNode>();
    root->right = right;
    EXPECT_CALL(*right, evaluate).Times(1);

    auto data = create_moldata(1, 1, 1, 1, 1);
    SelectionStack stack(root);
    stack.evaluate(*data, flags, 0);

    stack.push_flags(flags);
    SelectionFlags other = stack.pop_flags();
    ASSERT_THAT(flags.mask, NotNull());
    EXPECT_EQ(flags.mask, other.mask);
    ASSERT_THAT(flags.selected, NotNull());
    EXPECT_EQ(flags.selected, other.selected);
}

TEST(Selection, SelectionParser) {
    EXPECT_THROW(SelectionParser("not valid"), MolError);
    EXPECT_THROW(SEL_PARSER.parse("not valid"), MolError);
    EXPECT_TRUE(SEL_PARSER.parse("resid 1"));
}

TEST(Selection, BooleanParsing) {
    // Mock MolData
    auto data = create_moldata(5, 2, 1, 1, 1);
    std::shared_ptr<std::set<size_t>> selected;

    auto sel_tree = SEL_PARSER.parse("resid 1 or resid 2 or resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(2, 3, 4, 5, 6, 7));

    sel_tree = SEL_PARSER.parse("resid 1 and resid 2 and resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre());

    sel_tree = SEL_PARSER.parse("resid 1 or resid 2 and resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre());

    sel_tree = SEL_PARSER.parse("(resid 1 or resid 2) and resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre());

    sel_tree = SEL_PARSER.parse("resid 1 or (resid 2 and resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(2, 3));

    sel_tree = SEL_PARSER.parse("resid 1 and resid 2 or resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(6, 7));

    sel_tree = SEL_PARSER.parse("(resid 1 and resid 2) or resid 3");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(6, 7));

    sel_tree = SEL_PARSER.parse("resid 1 and (resid 2 or resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre());

    sel_tree = SEL_PARSER.parse("not resid 1 or resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(0, 1, 4, 5, 6, 7, 8, 9));

    sel_tree = SEL_PARSER.parse("resid 1 or not resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(0, 1, 2, 3, 6, 7, 8, 9));

    sel_tree = SEL_PARSER.parse("not resid 1 and resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left->left));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(4, 5));

    sel_tree = SEL_PARSER.parse("resid 1 and not resid 2");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(2, 3));

    sel_tree = SEL_PARSER.parse("resid 1 or not (resid 2 or resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(0, 1, 2, 3, 8, 9));

    sel_tree = SEL_PARSER.parse("resid 1 and not (resid 2 or resid 3)");
    ASSERT_TRUE(std::dynamic_pointer_cast<AndSelection>(sel_tree));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NotSelection>(sel_tree->right));
    ASSERT_TRUE(std::dynamic_pointer_cast<OrSelection>(sel_tree->right->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left->left));
    ASSERT_TRUE(std::dynamic_pointer_cast<NumPropSelection>(sel_tree->right->left->right));
    selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(2, 3));
}

TEST(Selection, NumPropParsing) {
    // Resid
    auto sel_tree = SEL_PARSER.parse("resid 0 or resid 2:4 or resid 6:7 9 11 or resid 30");
    // Nodes construction
    ASSERT_TRUE(sel_tree);
    ASSERT_TRUE(sel_tree->left);
    ASSERT_TRUE(sel_tree->right);
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left->left));
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->left->right));
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->left->right));
    EXPECT_TRUE(std::dynamic_pointer_cast<ResidSelection>(sel_tree->right));
    // Selection
    auto data = create_moldata(12, 1, 1, 1, 1);
    std::shared_ptr<std::set<size_t>> selected = evaluate_sel_tree(sel_tree, *data);
    EXPECT_THAT(*selected, ElementsAre(0, 2, 3, 4, 6, 7, 9, 11));
}

TEST(Selection, AtomSelector) {
    PDBFiles pdb;
    pdb.check();

    // Valid selections
    AtomSelector selector("resid 203:205", pdb.big);
    AtomSel water = selector.apply(0);
    EXPECT_THAT(water.indices(), ElementsAre(1807, 1808, 1809));
    EXPECT_EQ(water.frame(), 0);

    // No frame
    water = selector.apply(std::nullopt);
    EXPECT_THAT(water.indices(), ElementsAre(1807, 1808, 1809));
    EXPECT_EQ(water.frame(), 0);

    // Invalid selections
    selector = AtomSelector("resid 900:910", pdb.big);
    EXPECT_THAT(selector.apply(0).indices(), ElementsAre());

    // Errors
    EXPECT_THROW(selector.apply(1), MolError);
    EXPECT_THROW(AtomSelector("nonsense", pdb.big), MolError);
}
