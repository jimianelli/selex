-- Remove markdown source code blocks that appear when embedding child .qmd files.
-- Keep executable chunk code blocks (R/Python/etc.) intact.
function CodeBlock(el)
  if el.classes:includes("markdown") then
    return {}
  end
  return el
end
