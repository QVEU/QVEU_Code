# Markdown Cheat Sheet

Adapted from  [The Markdown Guide](https://www.markdownguide.org).

This Markdown cheat sheet provides a quick overview of all the Markdown syntax elements. It can’t cover every edge case, so if you need more information about any of these elements, refer to the reference guides for [basic syntax](https://www.markdownguide.org/basic-syntax) and [extended syntax](https://www.markdownguide.org/extended-syntax).

## Basic Syntax

These are the elements outlined in John Gruber’s original design document. All Markdown applications support these elements.

## Heading

\# H1
# H1

\## H2
## H2
\### H3
### H3

## Bold

```
**bold text**
```
**bold text**

## Italic

```
*italicized text*
```
*italicized text*

## Blockquote
```
> blockquote
```
> blockquote

## Ordered List
```
1. First item 
2. Second item 
3. Third item 
```

1. First item
2. Second item
3. Third item

## Unordered List

```
- First item
- Second item 
- Third item
```
- First item
- Second item
- Third item
---
## Code

```
`code`
```
`code`

## Horizontal Rule

```
---
```
---

### Link
```
[QVEU Website](http://qveu.bio)
```

[QVEU Website](http://qveu.bio)

### Image
```
![alt t-rex text](https://upload.wikimedia.org/wikipedia/commons/b/b2/Tyrannosaurus_skeletal_diagram.jpg)
```

![alt t-rex text](https://upload.wikimedia.org/wikipedia/commons/b/b2/Tyrannosaurus_skeletal_diagram.jpg)

# Extended Syntax

These elements extend the basic syntax by adding additional features. Not all Markdown applications support these elements.

## Table
```
| Syntax | Description |
| ----------- | ----------- |
| Header | Title |
| Paragraph | Text |
```

| Syntax | Description |
| ----------- | ----------- |
| Header | Title |
| Paragraph | Text |

## Fenced Code Block
```
'```
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```'
```
## Footnote

```
Here's a sentence with a footnote. [^1]
```
Here's a sentence with a footnote. [^1]

```
[^1]: This is the footnote. #(See below)
```
[^1]: This is the footnote.

## Heading ID

## Definition List
```
term
: definition
```
term
: definition

### Strikethrough

```
~~The world is flat.~~
```
~~The world is flat.~~

### Task List
```
- [x] Pass cells
- [x] Sequence Virus
- [ ] Save the World
```
- [x] Pass cells
- [x] Sequence Virus
- [ ] Save the World
### Emoji

```What an amazing discovery! :rocket:```
What an amazing discovery! :rocket:

(See also [Copying and Pasting Emoji](https://www.markdownguide.org/extended-syntax/#copying-and-pasting-emoji))

### Subscript

```
H<sub>2</sub>O
```
H<sub>2</sub>O

### Superscript
```
<sup>2</sup>
```
e=mc<sup>2</sup>
