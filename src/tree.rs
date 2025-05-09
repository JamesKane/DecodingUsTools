use crate::vendor::{HaplogroupTreeProvider, VariantProvider};

pub struct TreeLoader {
    provider: Box<dyn HaplogroupTreeProvider>,
    variant_parser: Box<dyn VariantProvider>,
}

impl TreeLoader {
    pub fn new(
        provider: Box<dyn HaplogroupTreeProvider>,
        variant_parser: Box<dyn VariantProvider>,
    ) -> Self {
        Self {
            provider,
            variant_parser,
        }
    }

    pub fn load_tree(&self, tree_type: TreeType) -> Result<HaplogroupTree, Box<dyn std::error::Error>> {
        self.provider.load_tree(tree_type)
    }
}